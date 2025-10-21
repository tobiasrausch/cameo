#ifndef CAMEO_H
#define CAMEO_H

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "edlib.h"
#include "util.h"
#include "pileup.h"

namespace cameo {

  struct CameoConfig {
    bool onlyCpG;
    bool combineStrands;
    bool hasBedFile;
    bool callPeaks;
    bool modifiedRegions;
    char peakStrand;
    char modCode;
    uint16_t minBaseQual;
    uint16_t minMapQual;
    uint32_t peakwin;
    uint32_t peakmedwin;
    uint32_t minCov;
    float minMod;
    float maxUnmod;
    float peakmethylth;
    boost::filesystem::path bedfile;
    boost::filesystem::path outfile;
    boost::filesystem::path bamfile;
    boost::filesystem::path genome;
    std::string sampleName;
  };
  
  template<typename TConfig>
  inline int32_t
  runCameo(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("cameo.prof");
#endif

    // Aggregate modifications
    screenMods(c);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    
    return 0;
  }

  int pileup(int argc, char **argv) {
    CameoConfig c;
   
    // Parameter
    std::string instag;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("base-qual,b", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(1), "min. base quality")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName), "sample name")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output pileup table")
     ;
    
    boost::program_options::options_description methyl("Methylation options");
    methyl.add_options()
      ("minmod,m", boost::program_options::value<float>(&c.minMod)->default_value(0.7), "min. mod probability threshold")
      ("maxunmod,u", boost::program_options::value<float>(&c.maxUnmod)->default_value(0.2), "max. mod probability threshold for unmodified")
      ("bedfile,f", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "report mods over input BED")
      ("peaks,k", "identify (de)methylated regions")
      ("cpg,p", "only CpG counts")
      ("combine,c", "combine strands")
      ;

    boost::program_options::options_description dmr("(De)methylated regions (requires -p)");
    dmr.add_options()
      ("movavg,a", boost::program_options::value<uint32_t>(&c.peakwin)->default_value(5), "rolling average window")
      ("medsize,e", boost::program_options::value<uint32_t>(&c.peakmedwin)->default_value(11), "rolling median window")
      ("strand,r", boost::program_options::value<char>(&c.peakStrand)->default_value('*'), "strand [+,-,*]")
      ("mincov,n", boost::program_options::value<uint32_t>(&c.minCov)->default_value(10), "min. coverage of each mod site")
      ("modletter,l", boost::program_options::value<char>(&c.modCode)->default_value('m'), "modification code [m,h]")
      ("thres,t", boost::program_options::value<float>(&c.peakmethylth)->default_value(0.8), "min. (un)modified fraction")
      ("modified,i", "identify peaks in modified signal instead of unmodified")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input file")
      ;
   
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(methyl).add(dmr).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(methyl).add(dmr);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: cameo " << argv[0] << " [OPTIONS] -g <ref.fa> <input.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }
    
    // Modification probability
    if (c.minMod < 0) c.minMod = 0;
    if (c.minMod > 1) c.minMod = 1;

    // CpGs
    if (vm.count("cpg")) c.onlyCpG = true;
    else c.onlyCpG = false;

    // Combine strands
    if (vm.count("combine")) c.combineStrands = true;
    else c.combineStrands = false;

    // Modified regions
    if (vm.count("modified")) c.modifiedRegions = true;
    else c.modifiedRegions = false;

    // Input BED file
    if (vm.count("bedfile")) c.hasBedFile = true;
    else c.hasBedFile = false;

    // Call DMRs
    if (vm.count("peaks")) {
      c.hasBedFile = true; // Peaks are stored as BED intervals
      c.callPeaks = true;
    } else c.callPeaks = false;
    
    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
    
    // Check input files
    if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
      std::cerr << "Alignment file is missing: " << c.bamfile.string() << std::endl;
      return 1;
    }
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
      return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
      return 1;
    }
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      std::string tname(hdr->target_name[refIndex]);
      if (!faidx_has_seq(fai, tname.c_str())) {
	std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	return 1;
      }
    }
    fai_destroy(fai);
    // Get sample name
    if (!vm.count("sample")) getSMTag(std::string(hdr->text), c.bamfile.stem().string(), c.sampleName);
    // Close files
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "cameo ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
   
    return runCameo(c);
  }

}

#endif
