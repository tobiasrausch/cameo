#ifndef PILEUP_H
#define PILEUP_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/format.hpp>

#include <htslib/sam.h>

#include "bed.h"
#include "util.h"

namespace cameo
{

  struct ModHit {
    int32_t pos;
    char code;
    uint8_t prob;
    bool rev;
    char base;

    ModHit(int32_t const p, char const c, uint8_t const r, bool const s, char const b) : pos(p), code(c), prob(r), rev(s), base(b) {}
  };

  template<typename TConfig>
  inline bool
  findPromoters(TConfig const& c) {
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Output table
    std::streambuf* buf;
    std::ofstream of;
    if(c.outfile.string() != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else buf = std::cout.rdbuf();
    std::ostream out(buf);
    out << "chr\tstart\tend\tsample\tmodbase\tstrand\tmodcount\tunmod\tothermod\tcoverage\tfail\tfracmod" << std::endl;

    // Load FASTA
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return false;
      } else fai = fai_load(c.genome.string().c_str());
    }
    
    // Keep track of avg. coverage
    typedef boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > TAccumulator;
    std::vector<TAccumulator> acc(c.files.size());
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Methylation scanning" << std::endl;

    // Parse BED
    typedef std::vector<BedEntry> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome scanRegions;
    if (c.hasBedFile) {
      int32_t nreg = _parseBedIntervals(c, scanRegions);
      if (nreg == 0) {
	std::cerr << "Error: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	return 1;
      } else {
	std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parsed " << nreg << " input regions." << std::endl;
      }
    }
    
    // Calling threshold
    uint8_t probModTh = (uint8_t) ((int) (c.minMod * 256));
    uint8_t probUnmodTh = (uint8_t) ((int) (c.maxUnmod * 256));
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Load chromosome
      std::string tname = hdr->target_name[refIndex];
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
      
      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {

	typedef uint16_t TCovValue;
	TCovValue maxval = std::numeric_limits<TCovValue>::max();
	std::vector<TCovValue> unmod_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> unmod_minus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> uncertain_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> uncertain_minus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> m_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> m_minus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> h_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> h_minus(hdr->target_len[refIndex], 0);
	
	// Read alignments
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  
	  // Keep secondary alignments
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  // Parse quality and seq
	  bool readRev = (bool) (rec->core.flag & BAM_FREVERSE);
	  typedef std::vector<uint8_t> TQuality;
	  TQuality quality;
	  quality.resize(rec->core.l_qseq);
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  uint8_t* qualptr = bam_get_qual(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	    quality[i] = qualptr[i];
	    sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  }
	  std::string fwdseq = sequence;
	  if (readRev) reverseComplement(fwdseq);

	  // Parse MM and ML tags (if present)
	  std::vector<ModHit> modhits;
	  uint8_t* mm_aux = bam_aux_get(rec, "MM");
	  if (mm_aux && *mm_aux == 'Z') {
	    const char* mmstr = (const char*)(mm_aux + 1);
	    std::string s(mmstr);
	    std::vector<std::string> tokens;
	    boost::split(tokens, s, boost::is_any_of(";"));
	    for(const auto& tok : tokens) {
	      if (tok.empty()) continue;
	      // parse <base><strand><modcode>,<pos>,<pos>,... e.g. C+m,0,5,...
	      std::size_t idx_tok = 0;
	      char base = tok[idx_tok++];
	      if (idx_tok >= tok.size()) continue;
	      char strand = tok[idx_tok++];
	      bool revMod = false;
	      if (strand == '-') revMod = true;
	      std::string mod_codes;
	      while ((idx_tok < tok.size()) && (tok[idx_tok] != ',')) {
		char ch = tok[idx_tok++];
		if (std::isalpha(static_cast<unsigned char>(ch))) mod_codes.push_back(ch); // ignore ?
	      }
	      // parse positions
	      if ((idx_tok < tok.size()) && (tok[idx_tok] == ',')) {
		std::string pos_str = tok.substr(idx_tok+1);
		if (!pos_str.empty()) {
		  std::vector<std::string> pos_tokens;
		  boost::split(pos_tokens, pos_str, boost::is_any_of(","));
		  int32_t current = -1;
		  for(const auto& ptoken : pos_tokens) {
		    if (ptoken.empty()) continue;
		    int32_t delta = std::stoi(ptoken);
		    current += (delta + 1);
		    for(char mc : mod_codes) modhits.push_back(ModHit(current, mc, 255, revMod, base));
		  }
		}
	      }
	    }
	  }

	  // Parse ML tag
	  uint8_t* ml_aux = bam_aux_get(rec, "ML");
	  if (ml_aux && *ml_aux == 'B') {
	    char subtype = *(char*)(ml_aux + 1);
	    if (subtype == 'C') {
	      const uint8_t* p = ml_aux + 2;
	      int32_t n = 0;
	      n = p[0] | (p[1]<<8) | (p[2]<<16) | (p[3]<<24);
	      const uint8_t* data = p + 4;
	      int32_t assign = std::min<int32_t>(n, (int32_t) modhits.size());
	      for(int32_t i = 0; i < assign; ++i) modhits[i].prob = data[i];
	    }
	  }
	  
	  // Find read position of As, Cs, ...
	  std::unordered_map<char, std::vector<int32_t>> basepos;
	  for (int32_t i = 0; i < (int32_t) fwdseq.size(); ++i) basepos[std::toupper(static_cast<unsigned char>(fwdseq[i]))].push_back(i);
	  
	  // Replace base positions with absolute read positions
	  std::unordered_map<int32_t, std::vector<ModHit> > modByPos;
	  for (const auto& mh : modhits) {
	    char ub = std::toupper(static_cast<unsigned char>(mh.base));
 	    char target_base = (mh.rev) ? complement_base(ub) : ub;
 	    auto it = basepos.find(target_base);
	    if (it == basepos.end()) continue;
	    const auto& occs = it->second;
	    if ((mh.pos < 0) || ((std::size_t)mh.pos >= occs.size())) continue;
	    // Absolute read pos is occs[mh.pos]
	    if ( (uint16_t) (quality[occs[mh.pos]]) >= c.minBaseQual) {
	      modByPos[occs[mh.pos]].push_back(ModHit(occs[mh.pos], mh.code, mh.prob, mh.rev, mh.base));
	    }
	  }
	  
	  // Parse cigar
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    uint32_t op = bam_cigar_op(cigar[i]);
	    uint32_t oplen = bam_cigar_oplen(cigar[i]);
	    if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	      for(uint32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
		if (rp >= hdr->target_len[refIndex]) break;
		// Filter for base quality
		if (quality[sp] >= c.minBaseQual) {
		  uint32_t lookUpPos = sp;
		  if (readRev) lookUpPos = rec->core.l_qseq - sp - 1;
		  auto it = modByPos.find(lookUpPos);
		  if (it != modByPos.end()) {
		    uint32_t unmod = 0;
		    uint32_t uncertain = 0;
		    std::set<char> modmatch;
		    for(const auto& mh : it->second) {
		      if (((mh.rev == readRev) && (sequence[sp] == mh.base)) || ((mh.rev != readRev) && (sequence[sp] == complement_base(mh.base)))) {
			if (mh.prob >= probModTh) modmatch.insert(mh.code);
			if (mh.prob > probUnmodTh) ++uncertain;
			else ++unmod;
		      }
		    }
		    if (modmatch.size() == 1) {
		      char code = *(modmatch.begin());
		      if ((code == 'm') || (code == 'M')) {
			if (readRev) {
			  if (m_minus[rp] < maxval) ++m_minus[rp];
			} else {
			  if (m_plus[rp] < maxval) ++m_plus[rp];
			}
		      } else if ((code == 'h') || (code == 'H')) {
			if (readRev) {
			  if (h_minus[rp] < maxval) ++h_minus[rp];
			} else {
			  if (h_plus[rp] < maxval) ++h_plus[rp];
			}
		      } else {
			std::cerr << "Warning: Unknown modification code! " << code << std::endl;
		      }
		    } else if ((modmatch.size() > 1) || (uncertain)) {
		      // Fail
		      if (readRev) {
			if (uncertain_minus[rp] < maxval) ++uncertain_minus[rp];
		      } else {
			if (uncertain_plus[rp] < maxval) ++uncertain_plus[rp];
		      }
		    } else {
		      if (readRev) {
			if (unmod_minus[rp] < maxval) ++unmod_minus[rp];
		      } else {
			if (unmod_plus[rp] < maxval) ++unmod_plus[rp];
		      }
		    }
		  }
		}
	      }
	    } else if (op == BAM_CDEL) {
	      rp += oplen;
	    } else if (op == BAM_CINS) {
	      sp += oplen;
	    } else if ((op == BAM_CSOFT_CLIP) || (op == BAM_CHARD_CLIP)) {
	      sp += oplen;
	    } else if (op == BAM_CREF_SKIP) {
	      rp += oplen;
	    } else {
	      std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
	
	// Summarize coverage for this chromosome
	for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  uint32_t total_cov = m_plus[pos] + h_plus[pos] + unmod_plus[pos] + m_minus[pos] + h_minus[pos] + unmod_minus[pos] + uncertain_plus[pos] + uncertain_minus[pos];
	  if (total_cov > 0) {
	    double failfrac = (double) (uncertain_plus[pos] + uncertain_minus[pos]) / (double) (total_cov);
	    acc[file_c](failfrac);
	  }
	}
	double avgfail = boost::accumulators::mean(acc[file_c]);
	std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << c.sampleName[file_c] << ", " << hdr->target_name[refIndex] << ", running failed fraction of mod calls=" << avgfail << std::endl;

	// Output percent modified per site
	if (c.hasBedFile) {
	  for(uint32_t i = 0; i < scanRegions[refIndex].size(); ++i) {
	    if ((scanRegions[refIndex][i].start >= 0) && (scanRegions[refIndex][i].end <= hdr->target_len[refIndex])) {
	      // Aggregate mod counts
	      uint32_t reg_unmod_plus = 0;
	      uint32_t reg_unmod_minus = 0;
	      uint32_t reg_uncertain_plus = 0;
	      uint32_t reg_uncertain_minus = 0;
	      uint32_t reg_m_plus = 0;
	      uint32_t reg_m_minus = 0;
	      uint32_t reg_h_plus = 0;
	      uint32_t reg_h_minus = 0;
	      for (uint32_t pos = scanRegions[refIndex][i].start; pos < scanRegions[refIndex][i].end; ++pos) {
		// CpG sites
		bool fwdCpG = true;
		bool revCpG = true;
		if (c.onlyCpG) {
		  fwdCpG = false;
		  revCpG = false;
		  if (pos + 1 < hdr->target_len[refIndex]) {
		    char b1 = std::toupper(seq[pos]);
		    char b2 = std::toupper(seq[pos+1]);
		    if ((b1 == 'C') && (b2 == 'G')) fwdCpG = true;
		  }
		  if (pos > 0) {
		    char b1 = std::toupper(seq[pos-1]);
		    char b2 = std::toupper(seq[pos]);
		    if ((b1 == 'C') && (b2 == 'G')) revCpG = true;
		  }
		}
		if (fwdCpG) {
		  reg_unmod_plus += unmod_plus[pos];
		  reg_uncertain_plus += uncertain_plus[pos];
		  reg_m_plus += m_plus[pos];
		  reg_h_plus += h_plus[pos];
		}
		if (revCpG) {
		  reg_unmod_minus += unmod_minus[pos];
		  reg_uncertain_minus += uncertain_minus[pos];				  
		  reg_m_minus += m_minus[pos];
		  reg_h_minus += h_minus[pos];
		}
	      }
	      // Unstranded
	      if (c.combineStrands) {
		uint32_t unstranded_cov = reg_m_plus + reg_h_plus + reg_unmod_plus + reg_m_minus + reg_h_minus + reg_unmod_minus;
		if (unstranded_cov > 0) {
		  double pct_h = double(reg_h_plus + reg_h_minus) / double(unstranded_cov);
		  double pct_m = double(reg_m_plus + reg_m_minus) / double(unstranded_cov);
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\th\t*\t" << (reg_h_plus + reg_h_minus) << "\t" << (reg_unmod_plus + reg_unmod_minus) << "\t" << (reg_m_plus + reg_m_minus) << "\t" << unstranded_cov << "\t" << (reg_uncertain_plus + reg_uncertain_minus) << "\t" << (boost::format("%1$.4f") % pct_h) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\tm\t*\t" << (reg_m_plus + reg_m_minus) << "\t" << (reg_unmod_plus + reg_unmod_minus) << "\t" << (reg_h_plus + reg_h_minus) << "\t" << unstranded_cov << "\t" << (reg_uncertain_plus + reg_uncertain_minus) << "\t" << (boost::format("%1$.4f") % pct_m) << std::endl;
		}
	      } else {
		// Plus strand
		uint32_t reg_cov_plus = reg_m_plus + reg_h_plus + reg_unmod_plus;
		if (reg_cov_plus) {
		  double pct_h_plus = double(reg_h_plus) / double(reg_cov_plus);
		  double pct_m_plus = double(reg_m_plus) / double(reg_cov_plus);
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\th\t+\t" << reg_h_plus << "\t" << reg_unmod_plus << "\t" << reg_m_plus << "\t" << reg_cov_plus << "\t" << reg_uncertain_plus << "\t" << (boost::format("%1$.4f") % pct_h_plus) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\tm\t+\t" << reg_m_plus << "\t" << reg_unmod_plus << "\t" << reg_h_plus << "\t" << reg_cov_plus << "\t" << reg_uncertain_plus << "\t" << (boost::format("%1$.4f") % pct_m_plus) << std::endl;
		}
		uint32_t reg_cov_minus = reg_m_minus + reg_h_minus + reg_unmod_minus;
		if (reg_cov_minus) {
		  double pct_h_minus = double(reg_h_minus) / double(reg_cov_minus);
		  double pct_m_minus = double(reg_m_minus) / double(reg_cov_minus);
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\th\t-\t" << reg_h_minus << "\t" << reg_unmod_minus << "\t" << reg_m_minus << "\t" << reg_cov_minus << "\t" << reg_uncertain_minus << "\t" << (boost::format("%1$.4f") % pct_h_minus) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << scanRegions[refIndex][i].start << "\t" << scanRegions[refIndex][i].end << "\t" << c.sampleName[file_c] << "\tm\t-\t" << reg_m_minus << "\t" << reg_unmod_minus << "\t" << reg_h_minus << "\t" << reg_cov_minus << "\t" << reg_uncertain_minus << "\t" << (boost::format("%1$.4f") % pct_m_minus) << std::endl;
		}
	      }
	    }
	  }
	} else {
	  // No BED input
	  for (uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	    // CpG sites
	    bool fwdCpG = true;
	    bool revCpG = true;
	    if (c.onlyCpG) {
	      fwdCpG = false;
	      revCpG = false;
	      if (pos + 1 < hdr->target_len[refIndex]) {
		char b1 = std::toupper(seq[pos]);
		char b2 = std::toupper(seq[pos+1]);
		if ((b1 == 'C') && (b2 == 'G')) fwdCpG = true;
	      }
	      if (pos > 0) {
		char b1 = std::toupper(seq[pos-1]);
		char b2 = std::toupper(seq[pos]);
		if ((b1 == 'C') && (b2 == 'G')) revCpG = true;
	      }
	    }

	    // Unstranded
	    if (c.combineStrands) {
	      if ((c.onlyCpG) && (fwdCpG) && (pos + 1 < hdr->target_len[refIndex])) {
		uint32_t unstranded_cov = m_plus[pos] + h_plus[pos] + unmod_plus[pos] + m_minus[pos+1] + h_minus[pos+1] + unmod_minus[pos+1];
		if (unstranded_cov > 0) {
		  double pct_h = double(h_plus[pos] + h_minus[pos+1]) / double(unstranded_cov);
		  double pct_m = double(m_plus[pos] + m_minus[pos+1]) / double(unstranded_cov);
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\th\t*\t" << (h_plus[pos] + h_minus[pos+1]) << "\t" << (unmod_plus[pos] + unmod_minus[pos+1]) << "\t" << (m_plus[pos] + m_minus[pos+1]) << "\t" << unstranded_cov << "\t" << (uncertain_plus[pos] + uncertain_minus[pos+1]) << "\t" << (boost::format("%1$.4f") % pct_h) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\tm\t*\t" << (m_plus[pos] + m_minus[pos+1]) << "\t" << (unmod_plus[pos] + unmod_minus[pos+1]) << "\t" << (h_plus[pos] + h_minus[pos+1]) << "\t" << unstranded_cov << "\t" << (uncertain_plus[pos] + uncertain_minus[pos+1]) << "\t" << (boost::format("%1$.4f") % pct_m) << std::endl;
		}
	      }
	    } else {
	      // Plus strand
	      if (fwdCpG) {
		uint32_t cov_plus = m_plus[pos] + h_plus[pos] + unmod_plus[pos];
		if (cov_plus > 0) {
		  double pct_h_plus = double(h_plus[pos]) / double(cov_plus);
		  double pct_m_plus = double(m_plus[pos]) / double(cov_plus);
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\th\t+\t" << h_plus[pos] << "\t" << unmod_plus[pos] << "\t" << m_plus[pos] << "\t" << cov_plus << "\t" << uncertain_plus[pos] << "\t" << (boost::format("%1$.4f") % pct_h_plus) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\tm\t+\t" << m_plus[pos] << "\t" << unmod_plus[pos] << "\t" << h_plus[pos] << "\t" << cov_plus << "\t" << uncertain_plus[pos] << "\t" << (boost::format("%1$.4f") % pct_m_plus) << std::endl;
		}
	      }
	      // Minus strand
	      if (revCpG) {
		uint32_t cov_minus = m_minus[pos] + h_minus[pos] + unmod_minus[pos];
		if (cov_minus > 0) {
		  double pct_h_minus = double(h_minus[pos]) / double(cov_minus);
		  double pct_m_minus = double(m_minus[pos]) / double(cov_minus);
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\th\t-\t" << h_minus[pos] << "\t" << unmod_minus[pos] << "\t" << m_minus[pos] << "\t" << cov_minus << "\t" << uncertain_minus[pos] << "\t" << (boost::format("%1$.4f") % pct_h_minus) << std::endl;
		  out << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\t" << c.sampleName[file_c] << "\tm\t-\t" << m_minus[pos] << "\t" << unmod_minus[pos] << "\t" << h_minus[pos] << "\t" << cov_minus << "\t" << uncertain_minus[pos] << "\t" << (boost::format("%1$.4f") % pct_m_minus) << std::endl;
		}
	      }
	    }
	  }
	}
      }
      // Free space
      if (seq != NULL) {
	free(seq);
	seq = NULL;
      }
    }

    // Close file
    if(c.outfile.string() != "-") of.close();

    // Clean-up
    if (fai) fai_destroy(fai);
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }

    return true;
  }
  
}

#endif
