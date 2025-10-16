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
	    auto it = base_occurrence_positions.find(ub);
	    if (it == base_occurrence_positions.end()) {
	      // no occurrences of that base in this read; skip
	      continue;
	    }
	    const auto& occs = it->second;
	    if (mh.pos < 0 || (std::size_t)mh.pos >= occs.size()) {
	      // occurrence index out of range -> skip
	      continue;
	    }
	    int32_t read_pos = occs[mh.pos]; // absolute read coordinate (sp)
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
	  bool readRev = (rec->core.flag & BAM_FREVERSE);
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
