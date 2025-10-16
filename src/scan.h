#ifndef SCAN_H
#define SCAN_H

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

#include "util.h"

#include <tuple>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cctype>

namespace cameo
{

  struct ModHit {
    int32_t pos;
    char code;
    uint8_t prob;
    char strand;
    char base;

    ModHit(int32_t const p, char const c, uint8_t const r, char const s, char const b) : pos(p), code(c), prob(r), strand(s), base(b) {}
  };


  static inline int32_t readpos_to_refpos_for_debug(bam1_t* rec, int32_t target_sp) {
    uint32_t rp = rec->core.pos; // 0-based reference position
    uint32_t sp = 0;
    uint32_t* cigar = bam_get_cigar(rec);
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
      uint32_t op = bam_cigar_op(cigar[i]);
      uint32_t oplen = bam_cigar_oplen(cigar[i]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
        for (uint32_t k = 0; k < oplen; ++k) {
          if (sp == (uint32_t)target_sp) return (int32_t)rp;
          ++rp; ++sp;
        }
      } else if (op == BAM_CDEL) {
        rp += oplen;
      } else if (op == BAM_CINS) {
        for (uint32_t k=0; k<oplen; ++k) {
          if (sp == (uint32_t)target_sp) {
            return (int32_t)rp;
          }
          ++sp;
        }
      } else if (op == BAM_CSOFT_CLIP) {
        for (uint32_t k=0; k<oplen; ++k) {
          if (sp == (uint32_t)target_sp) {
            return -1;
          }
          ++sp;
        }
      } else if (op == BAM_CHARD_CLIP) {
      } else if (op == BAM_CREF_SKIP) {
        rp += oplen;
      }
    }
    return -2;
  }

  // raw_modhits: vector of tuples (occurrence_index, modcode, prob, strand, base)
  static inline void debug_mm_for_read(bam1_t* rec, bam_hdr_t* hdr, const std::vector<std::tuple<int32_t,char,uint8_t,char,char>>& raw_modhits) {
    const char* qname = bam_get_qname(rec);
    bool readRev = (rec->core.flag & BAM_FREVERSE);
    uint8_t* mm_aux = bam_aux_get(rec, "MM");
    uint8_t* ml_aux = bam_aux_get(rec, "ML");
    std::string mm_raw = mm_aux ? std::string((char*)(mm_aux+1)) : std::string();
    std::string ml_repr = ml_aux ? "present" : "absent";

    std::cerr << "=== MM DEBUG ===\n";
    std::cerr << "QNAME: " << (qname?qname:"(null)") << " FLAG: " << rec->core.flag << " readRev: " << readRev << " POS: " << rec->core.pos << " CIGAR: ";
    // print cigar as simple string
    uint32_t* cigar = bam_get_cigar(rec);
    for (std::size_t i=0;i<rec->core.n_cigar;++i) {
      uint32_t op = bam_cigar_op(cigar[i]);
      uint32_t oplen = bam_cigar_oplen(cigar[i]);
      char oc = "MIDNSHP=XB"[op];
      std::cerr << oplen << oc;
    }
    std::cerr << "\n";
    std::cerr << "MM raw: " << mm_raw << " ML: " << ml_repr << "\n";

    // build read sequence string
    std::string seq;
    seq.resize(rec->core.l_qseq);
    uint8_t* seqptr = bam_get_seq(rec);
    for (int32_t i=0;i<rec->core.l_qseq;++i) seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr,i)];

    // For each raw modhit provide diagnostic mapping
    for (const auto& t : raw_modhits) {
      int32_t occ = std::get<0>(t);
      char code = std::get<1>(t);
      uint8_t prob = std::get<2>(t);
      char strand = std::get<3>(t);
      char base = std::get<4>(t);

      char ub = std::toupper(static_cast<unsigned char>(base));
      char target_base;
      if (strand == '+') target_base = ub;
      else {
        // complement (simple)
        switch (ub) {
        case 'A': target_base = 'T'; break;
        case 'T': target_base = 'A'; break;
        case 'C': target_base = 'G'; break;
        case 'G': target_base = 'C'; break;
        default: target_base = ub; break;
        }
      }

      std::vector<int32_t> occs;
      for (int32_t i=0;i<(int32_t)seq.size();++i)
        if (std::toupper(static_cast<unsigned char>(seq[i]))==target_base) occs.push_back(i);

      std::cerr << "ModHit token: base="<<base<<" strand="<<strand<<" code="<<code<<" occ_index(token)="<<occ<<" prob="<<(int)prob<<"\n";
      std::cerr << " target_base used to search read: " << target_base << " (occurrences in read="<<occs.size()<<")\n";
      if (occs.empty()) { std::cerr << "  -> No occurrences of target_base in read sequence\n"; continue; }

      // compute occ_index used by current code (mirroring for '-' tokens)
      std::size_t occ_index;
      if (strand == '+') occ_index = (std::size_t)occ;
      else occ_index = occs.size() - 1 - (std::size_t)occ;

      std::cerr << " occs list: ";
      for (auto x:occs) std::cerr << x << ",";
      std::cerr << "\n";
      std::cerr << " occ_index chosen = " << occ_index;
      if (occ_index < occs.size()) {
        std::cerr << " (read_pos/sp = " << occs[occ_index] << ")\n";
      } else {
        std::cerr << " (out-of-range occ_index)\n";
        continue;
      }
      int32_t sp = occs[occ_index];

      int32_t rp = readpos_to_refpos_for_debug(rec, sp);
      std::cerr << " mapped read_pos(sp="<<sp<<") -> ref_pos(rp) = ";
      if (rp >= 0) std::cerr << rp << " (0-based)\n"; else if (rp == -1) std::cerr << "soft-clipped (no ref mapping)\n"; else std::cerr << "out-of-range\n";

      // show how modOnRefPlus would be computed under both conventions
      bool modOnRefPlus_xor = ((strand == '+') != readRev);
      bool modOnRefPlus_direct = (strand == '+');
      std::cerr << " modOnRefPlus (XOR-with-readRev): " << modOnRefPlus_xor << "  (direct-token): " << modOnRefPlus_direct << "\n";
      std::cerr << "-----------------------\n";
    }
    std::cerr << "=== /MM DEBUG ===\n";
  }
  // --- /MM debug helper functions ---


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
    //std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Methylation scanning" << std::endl;

    bool onlyCpG = false;
    uint8_t probTh = (uint8_t) ((int) (0.5 * 256));
    
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
	std::vector<TCovValue> cov_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> cov_minus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> m_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> m_minus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> h_plus(hdr->target_len[refIndex], 0);
	std::vector<TCovValue> h_minus(hdr->target_len[refIndex], 0);
	
	// Read alignments
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  
	  // Keep secondary alignments
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  // Parse quality and seq
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
		    int32_t delta = 0;
		    delta = std::stoi(ptoken);
		    current += (delta + 1);
		    for(char mc : mod_codes) modhits.push_back(ModHit(current, mc, 255, strand, base));
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

	  // Debug invocation (controlled by env var CAMEO_DEBUG_MM)
	  {
	    const char* dbg_env = std::getenv("CAMEO_DEBUG_MM");
	    if (dbg_env) {
	      // parse optional limit from envvar, default 100
	      int limit = 100;
	      try { limit = std::stoi(std::string(dbg_env)); } catch(...) {}
	      static int debug_count = 0;
	      if (debug_count < limit) {
		// prepare raw vector of tuples (occurrence_index, code, prob, strand, base)
		std::vector<std::tuple<int32_t,char,uint8_t,char,char>> raw_modhits;
		raw_modhits.reserve(modhits.size());
		for (const auto& mh : modhits) raw_modhits.emplace_back(mh.pos, mh.code, mh.prob, mh.strand, mh.base);
		if (!raw_modhits.empty()) {
		  debug_mm_for_read(rec, hdr, raw_modhits);
		  ++debug_count;
		}
	      }
	    }
	  }

	  // Reverse read flag used later
	  bool readRev = (rec->core.flag & BAM_FREVERSE);
	  
	  // Reorder by read position
	  std::unordered_map<char, std::vector<int32_t>> base_occurrence_positions;
	  for (int32_t i = 0; i < (int32_t) sequence.size(); ++i) {
	    char b = std::toupper(static_cast<unsigned char>(sequence[i]));
	    base_occurrence_positions[b].push_back(i);
	  }
	  // Now create a new vector of ModHit with pos replaced by the absolute read position
	  std::vector<ModHit> adjusted_modhits;
	  adjusted_modhits.reserve(modhits.size());
	  for (const auto& mh : modhits) {
	    char ub = std::toupper(static_cast<unsigned char>(mh.base));
 	    char read_base = (mh.strand == '+') ? ub : complement_base(ub);
	    if (readRev) read_base = complement_base(read_base);
 	    auto it = base_occurrence_positions.find(read_base);
	    if (it == base_occurrence_positions.end()) continue;
	    const auto& occs = it->second;
	    if (mh.pos < 0 || (std::size_t)mh.pos >= occs.size()) continue;
	    std::size_t occ_index = (std::size_t) mh.pos;
 	    int32_t read_pos = occs[occ_index];
	    adjusted_modhits.push_back(ModHit(read_pos, mh.code, mh.prob, mh.strand, mh.base));
	  }
	  modhits.swap(adjusted_modhits);

	  // Build map from read position to modification hits
	  std::unordered_map<int32_t, std::vector<ModHit> > modByPos;
	  for(const auto& mh : modhits) modByPos[mh.pos].push_back(mh);
	  
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
		// Coverage
		if (readRev) {
		  if (cov_minus[rp] < maxval) ++cov_minus[rp];
		} else {
		  if (cov_plus[rp] < maxval) ++cov_plus[rp];
		}
		// Mods
		auto it = modByPos.find(sp);
		if (it != modByPos.end()) {
		  for(const auto& mh : it->second) {
		    if ((mh.prob == 255) || (mh.prob >= probTh)) {
		      bool modOnRefPlus = ((mh.strand == '+') != readRev);
		      if ((mh.code == 'm') || (mh.code == 'M')) {
			if (modOnRefPlus) {
			  if (m_plus[rp] < maxval) ++m_plus[rp];
			} else {
			  if (m_minus[rp] < maxval) ++m_minus[rp];
			}
		      } else if ((mh.code == 'h') || (mh.code == 'H')) {
			if (modOnRefPlus) {
			  if (h_plus[rp] < maxval) ++h_plus[rp];
			} else {
			  if (h_minus[rp] < maxval) ++h_minus[rp];
			}
		      } else {
			std::cerr << "Warning: Unknown modification code! " << mh.code << std::endl;
		      }
		    }
		  }
		}
	      }
	    } else if (op == BAM_CDEL) {
	      rp += oplen;
	    } else if (op == BAM_CINS) {
	      sp += oplen;
	    } else if (op == BAM_CSOFT_CLIP) {
	      sp += oplen;
	    } else if (op == BAM_CHARD_CLIP) {
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
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) acc[file_c](cov_plus[i] + cov_minus[i]);
	double sdcov = sqrt(boost::accumulators::variance(acc[file_c]));
	double avgcov = boost::accumulators::mean(acc[file_c]);
	//std::cerr << c.files[file_c].string() << ", " << hdr->target_name[refIndex] << ", meancov=" << avgcov << ", sdcov=" << sdcov << std::endl;

	// Output percent modified per site
	// Header: file,chrom,1-based-pos,refbase,cov_plus,cov_minus,mod_m_plus,mod_m_minus,mod_h_plus,mod_h_minus,percent_modified(+),percent_modified(-),percent_modified(collapsed)
	for (uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  // optionally count only CpG sites
	  if (onlyCpG) {
	    if (pos + 1 >= hdr->target_len[refIndex]) continue;
	    char b1 = std::toupper(seq[pos]);
	    char b2 = std::toupper(seq[pos+1]);
	    if (!(b1 == 'C' && b2 == 'G')) continue;
	  }
	  int32_t total_plus = h_plus[pos] + m_plus[pos];
	  int32_t total_minus = h_minus[pos] + m_minus[pos];
	  double pct_h_minus = ((cov_minus[pos] > 0) ? 100.0 * double(h_minus[pos]) / double(cov_minus[pos]) : 0.0);
	  double pct_h_plus = ((cov_plus[pos] > 0) ? 100.0 * double(h_plus[pos]) / double(cov_plus[pos]) : 0.0);
	  double pct_m_minus = ((cov_minus[pos] > 0) ? 100.0 * double(m_minus[pos]) / double(cov_minus[pos]) : 0.0);
	  double pct_m_plus = ((cov_plus[pos] > 0) ? 100.0 * double(m_plus[pos]) / double(cov_plus[pos]) : 0.0);
	  if (total_minus) std::cerr << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\th\t" << cov_minus[pos] << "\t-\t" << pos << "\t" << (pos + 1) << "\t255,0,0\t" << cov_minus[pos] << "\t" << (boost::format("%1$.2f") % pct_h_minus) << "\t" << h_minus[pos] << std::endl;
	  if (total_plus) std::cerr << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\th\t" << cov_plus[pos] << "\t+\t" << pos << "\t" << (pos + 1) << "\t255,0,0\t" << cov_plus[pos] << "\t" << (boost::format("%1$.2f") % pct_h_plus) << "\t" << h_plus[pos] << std::endl;
	  if (total_minus) std::cerr << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\tm\t" << cov_minus[pos] << "\t-\t" << pos << "\t" << (pos + 1) << "\t255,0,0\t" << cov_minus[pos] << "\t" << (boost::format("%1$.2f") % pct_m_minus) << "\t" << m_minus[pos] << std::endl;
	  if (total_plus) std::cerr << hdr->target_name[refIndex] << "\t" << pos << "\t" << (pos + 1) << "\tm\t" << cov_plus[pos] << "\t+\t" << pos << "\t" << (pos + 1) << "\t255,0,0\t" << cov_plus[pos] << "\t" << (boost::format("%1$.2f") % pct_m_plus) << "\t" << m_plus[pos] << std::endl;
	}
      }
      // Free space
      if (seq != NULL) {
	free(seq);
	seq = NULL;
      }
      exit(-1);
    }

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
