#ifndef MOTIF_H
#define MOTIF_H

#include <limits>

#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"

namespace cameo
{

  struct MotifConfig {
    int32_t nchr;
    boost::filesystem::path outfile;
    boost::filesystem::path motiffile;
    boost::filesystem::path vcffile;
    boost::filesystem::path bamfile;
    boost::filesystem::path genome;
    std::string sampleName;
  };


  struct Motif {
    typedef boost::multi_array<uint16_t, 2> T2DArray;
    T2DArray matrix;

    std::string matrixId;
    std::string symbol;
  };

  template<typename TConfig>
  inline bool
  parseJaspar(TConfig const& c, std::vector<Motif>& pfms) {
    // Parse JASPAR
    std::ifstream file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.motiffile)) {
      file.open(c.motiffile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else file.open(c.motiffile.string().c_str(), std::ios_base::in);
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    int32_t acgt = 0;
    int32_t id = 0;
    while(std::getline(instream, gline)) {
      // Header
      if ((gline.size()) && (gline[0] == '>')) {
	id = pfms.size();
	pfms.resize(id+1);
	gline = gline.substr(1);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  pfms[id].matrixId = *tokIter++;
	  pfms[id].symbol = "NA";
	  if (tokIter != tokens.end()) {
	    pfms[id].symbol = *tokIter++;
	  }
	}
	acgt = 0;
      } else {
	if ((gline.size()) && ((gline[0] == 'A') || (gline[0] == 'C') || (gline[0] == 'G') || (gline[0] == 'T'))) {
	  // JASPAR format
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep("[");
	  Tokenizer tokens(gline, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if ((tokIter!=tokens.end()) && (++tokIter != tokens.end())) {
	    gline = *tokIter;
	    boost::char_separator<char> sep2("]");
	    Tokenizer tokens2(gline, sep2);
	    Tokenizer::iterator tokIter2 = tokens2.begin();
	    if (tokIter2 != tokens2.end()) {
	      gline = *tokIter2;
	    } else {
	      std::cerr << "JASPAR cannot be parsed!" << std::endl;
	      return false;
	    }
	  } else {
	    std::cerr << "JASPAR cannot be parsed!" << std::endl;
	    return false;
	  }
	}

	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	if (acgt == 0) { 
	  int32_t lenMotif = 0;
	  for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter) ++lenMotif;
	  pfms[id].matrix.resize(boost::extents[4][lenMotif]);
	}
	uint32_t col = 0;
	for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter, ++col) pfms[id].matrix[acgt][col] = boost::lexical_cast<int32_t>(*tokIter);

	// Debug code
	if (acgt == 3) {
	  std::cout << ">" << pfms[id].matrixId << ',' << pfms[id].symbol << std::endl;
	  for(uint32_t i = 0; i < pfms[id].matrix.shape()[0]; ++i) {
	    for(uint32_t j = 0; j < pfms[id].matrix.shape()[1]; ++j) {
	      std::cerr << pfms[id].matrix[i][j] << ',';
	    }
	    std::cerr << std::endl;
	  }
	}
	
	++acgt;
      }

    }
    // Clean-up
    dataIn.pop();
    if (is_gz(c.motiffile)) dataIn.pop();
    return true;
  }

  template<typename TConfig>
  inline int
  runMotif(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("cameo.prof");
#endif
    
    // Parse somatic variants
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse VCF/BCF variants" << std::endl;
    std::vector<std::vector<Variant>> variants_by_chr(c.nchr);
    parseVcfVariants(c, variants_by_chr);

    // Parse motifs
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse motifs" << std::endl;
    std::vector<Motif> mtf;
    parseJaspar(c, mtf);

#ifdef PROFILE
    ProfilerStop();
#endif

    // Done
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    
    return 0;
  }


  int motif(int argc, char **argv) {
    MotifConfig c;
   
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("vcffile,f", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "VCF/BCF file with somatic variants")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName), "sample name")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output file")
     ;
    
    boost::program_options::options_description phasing("Mofif options");
    phasing.add_options()
      ("motif,m", boost::program_options::value<boost::filesystem::path>(&c.motiffile), "mofif file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input file")
      ;
   
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(phasing).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(phasing);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("vcffile"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: cameo " << argv[0] << " [OPTIONS] -g <ref.fa> -f <input.vcf> -m <motif.txt> <input.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }

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
    c.nchr = hdr->n_targets;
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

    // Check input VCF/BCF file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);
    if (!ibcffile) {
      std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
      return 1;
    }
    if (!bcfhdr) {
      std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(bcfhdr);
    bcf_close(ibcffile);
    
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
   
    return runMotif(c);
  }
  
}

#endif
