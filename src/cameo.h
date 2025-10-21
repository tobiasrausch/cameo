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
    uint16_t minBaseQual;
    uint16_t minMapQual;
    int32_t nchr;
    float minMod;
    float maxUnmod;
    boost::filesystem::path bedfile;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    std::vector<std::string> sampleName;
  };
  
  template<typename TConfig>
  inline int32_t
  runCameo(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("cameo.prof");
#endif

    // Search promoters
    findPromoters(c);
    
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
      ("base-qual,a", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(1), "min. base quality")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output pileup table")
     ;
    
    boost::program_options::options_description methyl("Methylation options");
    methyl.add_options()
      ("minmod,m", boost::program_options::value<float>(&c.minMod)->default_value(0.7), "min. mod probability threshold")
      ("maxunmod,u", boost::program_options::value<float>(&c.maxUnmod)->default_value(0.1), "max. mod probability threshold for unmodified")
      ("bedfile,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "report mods over input BED")
      ("cpg,p", "only CpG counts")
      ("combine,c", "combine strands")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
      ;
   
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(methyl).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(methyl);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: cameo " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
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

    // Input BED file
    if (vm.count("bedfile")) c.hasBedFile = true;
    else c.hasBedFile = false;
    
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
    c.sampleName.resize(c.files.size());
    c.nchr = 0;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      if (!c.nchr) c.nchr = hdr->n_targets;
      else {
	if (c.nchr != hdr->n_targets) {
	  std::cerr << "BAM files have different number of chromosomes!" << std::endl;
	  return 1;
	}
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
      std::string sampleName = "unknown";
      getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName);
      c.sampleName[file_c] = sampleName;
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
    checkSampleNames(c);

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
