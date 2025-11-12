#ifndef PHASE_H
#define PHASE_H

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "edlib.h"

namespace cameo {
  
  struct PhaseConfig {
    bool hasTaggedBam;
    uint16_t minBaseQual;
    uint16_t minMapQual;
    uint32_t maxThreads;
    uint32_t maxReadsPerChunk;
    uint32_t maxIter;
    uint32_t minReadsPerBlock;
    int32_t nchr;
    int32_t chunkSize;
    int32_t chunkOverlap;
    boost::filesystem::path outfile;
    boost::filesystem::path taggedfile;
    boost::filesystem::path vcffile;
    boost::filesystem::path bamfile;
    boost::filesystem::path genome;
    std::string sampleName;
  };

  struct Variant {
    int32_t pos;
    int32_t ps;
    int32_t hap;
    std::string ref;
    std::string alt;

    Variant() : pos(0), ps(0), hap(0), ref(""), alt("") {}
    Variant(int32_t const p, std::string const& r, std::string const& a) : pos(p), ps(0), hap(0), ref(r), alt(a) {}

    bool operator<(const Variant& v2) const {
      return (pos<v2.pos);
    }
  };

  struct ReadInfo {
    boost::dynamic_bitset<> mask;
    boost::dynamic_bitset<> alt;
    std::vector<uint8_t> bq;
    
    ReadInfo() : mask(), alt(), bq() {}
    explicit ReadInfo(int const numVars) {
      mask.resize(numVars, false);
      alt.resize(numVars, false);
      bq.resize(numVars, 0);
    }    
  };

  struct DPPath {
    double score;
    int hp; // 0 for H1, 1 for H2

    DPPath() : score(0.0), hp(0) {}
  };

  inline double
  _phredToErr(int q, int cap=60) {
    if (q < 0) q = 0;
    if (q > cap) q = cap;
    return std::pow(10.0, -q / 10.0);
  }
  
  inline int
  _bestAllele(std::string const& sequence, int sp, char const* seq, int gp, std::string const& REF, std::string const& ALT) {
    int buffer = 50;
    int maxEdit = 10;
    int32_t diff = std::abs((int) ALT.size() - (int) REF.size());
    if ((sp < buffer) || (gp < buffer) || (sp + diff + buffer > (int) sequence.size())) return -1;
    std::string ALTEXT;
    std::string REFEXT;
    if (REF.size() <= ALT.size()) {
      // Insertion
      ALTEXT = ALT + std::string(seq + gp + REF.size() - buffer, seq + gp + REF.size() + buffer);
      REFEXT = REF + std::string(seq + gp + REF.size() - buffer, seq + gp + REF.size() + diff + buffer);
    } else {
      // Deletion
      ALTEXT = ALT + std::string(seq + gp + REF.size() - buffer, seq + gp + REF.size() + diff + buffer);
      REFEXT = REF + std::string(seq + gp + REF.size() - buffer, seq + gp + REF.size() + buffer);
    }
    std::string query = sequence.substr(sp - buffer, diff + 2*buffer + 1);
    EdlibAlignResult rRef = edlibAlign(REFEXT.c_str(), REFEXT.size(), query.c_str(), query.size(), edlibNewAlignConfig(maxEdit, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    EdlibAlignResult rAlt = edlibAlign(ALTEXT.c_str(), ALTEXT.size(), query.c_str(), query.size(), edlibNewAlignConfig(maxEdit, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    int dRef = rRef.editDistance;
    int dAlt = rAlt.editDistance;
    edlibFreeAlignResult(rRef);
    edlibFreeAlignResult(rAlt);
    //std::cerr << ALTEXT.size() << ',' << REFEXT.size() << ',' << query.size() << ',' << dRef << ',' << dAlt << std::endl;
    
    // Confident ALT or REF call?
    if ((dRef < 0) && (dAlt < 0)) return -1;
    if ((dRef >= 0) && (dRef <= maxEdit) && ((dAlt < 0) || (dRef < dAlt))) return 0;
    if ((dAlt >= 0) && (dAlt <= maxEdit) && ((dRef < 0) || (dAlt < dRef))) return 1;
    return -1;
  }

  template<typename TConfig>
  inline void
  parseVcfVariants(TConfig const& c, std::vector<std::vector<Variant>>& variants_by_chr) {
    // Load BAM file
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load BCF file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

    // Fetch sample
    int32_t sampleIndex = -1;
    for (int32_t i = 0; i < bcf_hdr_nsamples(bcfhdr); ++i) {
      if (bcfhdr->samples[i] == c.sampleName) sampleIndex = i;
    }
    if (sampleIndex < 0) {
      if (bcf_hdr_nsamples(bcfhdr) > 0) {
	std::cerr << "Warning: None of the samples in the VCF/BCF file matches " << c.sampleName << "! Taking first sample in VCF/BCF for phasing: " << bcfhdr->samples[0] << std::endl;
	sampleIndex = 0;
      } else {
	std::cerr << "Warning: No genotypes are present in the VCF/BCF file. All variants are considered as het. bi-allelic variants." << std::endl;
      }	
    }

    // Parse BCF
    faidx_t* fai = fai_load(c.genome.string().c_str());
    char* seq = NULL;
    bcf1_t* rec = bcf_init();
    int32_t kept = 0;
    int32_t total = 0;
    int32_t lastRefIdx = -1;
    int32_t lastpos = -1;
    int32_t refIndex = -1;
    while (bcf_read(ibcffile, bcfhdr, rec) == 0) {
      ++total;
      bcf_unpack(rec, BCF_UN_STR);
      
      // New chromosome?
      if (rec->rid != lastRefIdx) {
	// Match rid to BAM tid
	std::string chrom = bcf_hdr_id2name(bcfhdr, rec->rid);
	refIndex = bam_name2id(hdr, chrom.c_str());
	if (refIndex < 0) {
	  std::cerr << "Warning: Chromosome does not exist in BAM file: " << chrom << std::endl;
	  continue;
	}
	if (seq != NULL) {
	  free(seq);
	  seq = NULL;
	}
	int32_t seqlen = -1;
	seq = faidx_fetch_seq(fai, chrom.c_str(), 0, hdr->target_len[refIndex], &seqlen);
	if (seqlen <= 0) std::cerr << "Warning: Chromosome does not exist in FASTA genome file: " << chrom << std::endl;
	lastpos = -1;
	lastRefIdx = rec->rid;
      }
      
      if (rec->n_allele == 2) {
	if (!rec->d.allele[0] || !rec->d.allele[1]) continue;
	int32_t ngt_arr=0;
	int32_t *gt_arr=NULL;
	if (bcf_get_genotypes(bcfhdr, rec, &gt_arr, &ngt_arr) < 0) { if (gt_arr) free(gt_arr); continue; }
	int a0 = bcf_gt_is_missing(gt_arr[sampleIndex*2]) ? -1 : bcf_gt_allele(gt_arr[sampleIndex*2]);
	int a1 = bcf_gt_is_missing(gt_arr[sampleIndex*2 + 1]) ? -1 : bcf_gt_allele(gt_arr[sampleIndex*2 + 1]);
	free(gt_arr);
	if ((a0 == a1) || (a0 < 0) || (a1 < 0) || (!((a0==0 && a1==1) || (a0==1 && a1==0)))) continue;
	if (rec->pos != lastpos) {
	  std::string rAllele = std::string(rec->d.allele[0]);
	  std::string aAllele = std::string(rec->d.allele[1]);
	  if (rAllele != boost::to_upper_copy(std::string(seq + rec->pos, seq + rec->pos + rAllele.size()))) {
	    std::cerr << "Warning: REF allele mismatch at " << hdr->target_name[refIndex] << ":" << (rec->pos + 1) << std::endl;
	  } else {
	    if (aAllele[0] != '*') {
	      variants_by_chr[refIndex].push_back(Variant(rec->pos, rAllele, aAllele));
	      ++kept;
	    }
	  }
	  lastpos = rec->pos;
	}
      }
    }
    bcf_destroy(rec);
    if (seq != NULL) free(seq);
    
    // Clean-up
    fai_destroy(fai);
    bcf_hdr_destroy(bcfhdr);
    bcf_close(ibcffile);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    // Sort variants by position
    for(uint32_t rIdx = 0; rIdx < variants_by_chr.size(); ++rIdx) {
      std::sort(variants_by_chr[rIdx].begin(), variants_by_chr[rIdx].end());
    }

    // Done
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parsed " << total << " variants, kept " << kept << " het. bi-allelic variants." << std::endl;
  }
  
  template<typename TConfig, typename TPosToLocalMap>
  inline void
  collectReadsPerChunk(TConfig const& c, int refIndex, std::vector<Variant> const&chunk_vars, TPosToLocalMap& pos_to_local, std::unordered_map<std::size_t, ReadInfo>& out_reads, int chunk_genomic_start, int chunk_genomic_end) {
    typedef std::unordered_map<std::size_t, ReadInfo> TReadMap;
    TReadMap readTmp;
    out_reads.clear();

    // Any variants to phase?
    if (!chunk_vars.empty()) {
      samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
      hts_set_fai_filename(samfile, c.genome.string().c_str());
      hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
      bam_hdr_t* hdr = sam_hdr_read(samfile);

      // Load reference sequence
      faidx_t* fai = fai_load(c.genome.string().c_str());
      std::string chrName(hdr->target_name[refIndex]);
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // Parse BAM
      hts_itr_t* itr = sam_itr_queryi(idx, refIndex, chunk_genomic_start, chunk_genomic_end);
      bam1_t* rec = bam_init1();
      int numVars = chunk_vars.size();
      while (sam_itr_next(samfile, itr, rec) >= 0) {
        if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
        if ((rec->core.l_qseq <= 0) || (rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t const* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	// Get base qualities
	typedef std::vector<uint8_t> TQuality;
	TQuality quality;
	quality.resize(rec->core.l_qseq);
	uint8_t* qualptr = bam_get_qual(rec);
	for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
	
	// Add ReadInfo
	std::size_t seed = hash_lr(bam_get_qname(rec));
	if (readTmp.find(seed) == readTmp.end()) readTmp[seed] = ReadInfo(numVars);
	ReadInfo &ri = readTmp[seed];

	// Parse cigar
        int gp = rec->core.pos;
        int sp = 0;
	uint32_t const* cigar = bam_get_cigar(rec);
        for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
	  if ((bam_cigar_op(cigar[ci]) == BAM_CMATCH) || (bam_cigar_op(cigar[ci]) == BAM_CEQUAL) || (bam_cigar_op(cigar[ci]) == BAM_CDIFF)) {
	    for (uint32_t i = 0; i < bam_cigar_oplen(cigar[ci]); ++i, ++gp, ++sp) {
	      if ((gp < chunk_genomic_start) || (gp >= chunk_genomic_end)) continue;
	      auto itpos = pos_to_local.find(gp);
	      if (itpos == pos_to_local.end()) continue;
	      int local = itpos->second;
	      if ((sp < 0) || (sp >= rec->core.l_qseq)) continue;
	      if (quality[sp] < c.minBaseQual) continue;
	      std::string const& REF = chunk_vars[local].ref;
	      std::string const& ALT = chunk_vars[local].alt;
	      int call = -1;
	      if ((ALT.size() == 1) && (REF.size() == 1)) {
		// SNP
		if ((sequence[sp] == REF[0]) && (sequence[sp] != ALT[0])) call = 0;
		else if ((sequence[sp] != REF[0]) && (sequence[sp] == ALT[0])) call = 1;
	      } else {
		if (REF.size() <= ALT.size()) {
		  // Insertion or MNP
		  int32_t diff = ALT.size() - REF.size();
		  std::string REFEXT = REF + std::string(seq + gp + REF.size(), seq + gp + REF.size() + diff);
		  if ((sequence.substr(sp, ALT.size()) != ALT) && (sequence.substr(sp, REFEXT.size()) == REFEXT)) call = 0;
		  else if ((sequence.substr(sp, ALT.size()) == ALT) && (sequence.substr(sp, REFEXT.size()) != REFEXT)) call = 1;
		  else call = _bestAllele(sequence, sp, seq, gp, REF, ALT);
		} else {
		  // Deletion
		  int32_t diff = REF.size() - ALT.size();
		  std::string ALTEXT = ALT + std::string(seq + gp + REF.size(), seq + gp + REF.size() + diff);
		  if ((sequence.substr(sp, ALTEXT.size()) != ALTEXT) && (sequence.substr(sp, REF.size()) == REF)) call = 0;
		  else if ((sequence.substr(sp, ALTEXT.size()) == ALTEXT) && (sequence.substr(sp, REF.size()) != REF)) call = 1;
		  else call = _bestAllele(sequence, sp, seq, gp, REF, ALT);
		}
	      }
	      if (call != -1) {
		ri.mask[local] = true;
		ri.alt[local] = call;
		ri.bq[local] = quality[sp];
	      }
	    }
	  } else if (bam_cigar_op(cigar[ci]) == BAM_CINS) {
	    sp += bam_cigar_oplen(cigar[ci]);
	  } else if ((bam_cigar_op(cigar[ci]) == BAM_CDEL) || (bam_cigar_op(cigar[ci]) == BAM_CREF_SKIP)) {
	    gp += bam_cigar_oplen(cigar[ci]);
	  } else if (bam_cigar_op(cigar[ci]) == BAM_CSOFT_CLIP) {
	    sp += bam_cigar_oplen(cigar[ci]);
	  } else if (bam_cigar_op(cigar[ci]) == BAM_CHARD_CLIP) {
	    // No shift for sp
	  }
        }
      }
      if (seq != NULL) free(seq);
      fai_destroy(fai);
      bam_destroy1(rec);
      hts_itr_destroy(itr);
      hts_idx_destroy(idx);
      bam_hdr_destroy(hdr);
      sam_close(samfile);
    
      // Remove reads that do not overlap >=2 variants and downsample (if needed)
      std::vector<std::pair<uint32_t, std::size_t> > varCounts;
      for(typename TReadMap::iterator it = readTmp.begin(); it != readTmp.end(); ++it) varCounts.push_back(std::make_pair(it->second.mask.count(), it->first));
      std::sort(varCounts.begin(), varCounts.end(), std::greater<>());
      uint32_t rcount = 0;
      for(uint32_t i = 0; i < varCounts.size(); ++i) {
	if (varCounts[i].first >= 2) {
	  if (rcount < c.maxReadsPerChunk) {
	    out_reads[varCounts[i].second] = readTmp[varCounts[i].second];
	    ++rcount;
	  }
	}
      }
    }
  }

  template<typename TConfig>
  inline void
  solveChunk(TConfig const& c, std::vector<Variant> const& chunk_vars, std::unordered_map<std::size_t, ReadInfo> const& reads_map, std::vector<int>& out_hap, std::vector<int>& phase_block_id) {
    if (chunk_vars.empty()) return;
    int B = (int)chunk_vars.size();
    out_hap.assign(B, -1);
    phase_block_id.assign(B, 0);

    // Map reads to integer
    std::vector<std::size_t> read_hashes;
    for (auto const& kv : reads_map) read_hashes.push_back(kv.first);
    int R = read_hashes.size();
    if ((R < (int) c.minReadsPerBlock) || (B <= 1)) return;

    // Split variants into blocks
    std::vector<int> spanCount(B - 1, 0);
    for (int r = 0; r < R; ++r) {
      ReadInfo const& ri = reads_map.at(read_hashes[r]);
      int last = -1;
      for (int i = 0; i < B; ++i) {
        if (ri.mask[i]) {
          if (last != -1) {
            for (int j = last; j < i; ++j) ++spanCount[j];
          }
          last = i;
        }
      }
    }
    std::vector<std::pair<int, int>> blocks;
    int block_start = 0;
    for (int i = 0; i < B - 1; ++i) {
      if (spanCount[i] < (int) c.minReadsPerBlock) {
        blocks.push_back(std::make_pair(block_start, i));
        block_start = i + 1;
      }
    }
    blocks.push_back(std::make_pair(block_start, B - 1));

    // Process blocks
    int current_ps_id = chunk_vars[0].pos;
    if (current_ps_id == 0) current_ps_id = 1;
    for (auto const& blk : blocks) {
      int s = blk.first;
      int e = blk.second;
      int block_size = e - s + 1;
      if (block_size < 2) continue;

      // Collect reads
      std::vector<int> block_read_indices;
      for (int r = 0; r < R; ++r) {
        ReadInfo const& ri = reads_map.at(read_hashes[r]);
        for (int i = s; i <= e; ++i) {
          if (ri.mask[i]) {
            block_read_indices.push_back(r);
            break;
          }
        }
      }
      if (block_read_indices.size() < c.minReadsPerBlock) continue;

      // Init HP1 (partition 0) and HP2 (partition 1)
      std::vector<int> partition(block_read_indices.size(), -1);
      partition[0] = 0;
      for (uint32_t i = 1; i < block_read_indices.size(); ++i) {
        ReadInfo const& ri = reads_map.at(read_hashes[block_read_indices[i]]);
        double score0 = 0;
	double score1 = 0;
        for (uint32_t j = 0; j < i; ++j) {
	  if (partition[j] == -1) continue;
	  ReadInfo const& rj = reads_map.at(read_hashes[block_read_indices[j]]);
	  double agree = 0;
	  double disagree = 0;
	  for (int k = s; k <= e; ++k) {
	    if (ri.mask[k] && rj.mask[k]) {
	      if (ri.alt[k] == rj.alt[k]) ++agree;
	      else ++disagree;
	    }
	  }
	  if (agree + disagree > 0) {
	    if (partition[j] == 0) score0 += disagree / (agree + disagree);
	    else score1 += disagree / (agree + disagree);
	  }
        }
        partition[i] = (score0 <= score1) ? 0 : 1;
      }

      // Iterative refinement
      std::vector<int> best_partition(block_read_indices.size(), 0);
      double min_mec_score = std::numeric_limits<double>::max();
      for (int iter = 0; iter < (int)c.maxIter; ++iter) {
        // Compute consensus
        std::vector<int> h1(block_size, -1);
	std::vector<int> h2(block_size, -1);
        for (uint32_t i = 0; i < (uint32_t) block_size; ++i) {
	  // H1
          int ref_c = 0;
	  int alt_c = 0;
          for (uint32_t r_idx = 0; r_idx < block_read_indices.size(); ++r_idx) {
            if (partition[r_idx] == 0) {
              ReadInfo const& ri = reads_map.at(read_hashes[block_read_indices[r_idx]]);
              if (ri.mask[s + i]) {
                if (ri.alt[s + i]) ++alt_c;
		else ++ref_c;
              }
            }
          }
          if (ref_c + alt_c > 0) h1[i] = (alt_c > ref_c) ? 1 : 0;
	  // H2
          ref_c = 0;
	  alt_c = 0;
          for (uint32_t r_idx = 0; r_idx < block_read_indices.size(); ++r_idx) {
            if (partition[r_idx] == 1) {
              ReadInfo const& ri = reads_map.at(read_hashes[block_read_indices[r_idx]]);
              if (ri.mask[s + i]) {
                if (ri.alt[s + i]) ++alt_c;
		else ++ref_c;
              }
            }
          }
          if (ref_c + alt_c > 0) h2[i] = (alt_c > ref_c) ? 1 : 0;
	}

        // DP
        std::vector<int> new_partition(block_read_indices.size());
        double current_mec_score = 0;
        bool changed = false;
        for (uint32_t r_idx = 0; r_idx < block_read_indices.size(); ++r_idx) {
          ReadInfo const& ri = reads_map.at(read_hashes[block_read_indices[r_idx]]);
          std::vector<DPPath> dp(block_size, DPPath());
          double switch_penalty = std::log(0.001);
          for (int i = 0; i < block_size; ++i) {
	    if (!ri.mask[s + i]) {
	      if (i > 0) dp[i] = dp[i-1];
	      continue;
	    }
	    int allele = ri.alt[s+i];
	    double err_prob = _phredToErr(ri.bq[s+i]);
	    double cost1 = ((h1[i] == -1) || (h1[i] == allele)) ? err_prob : (1.0 - err_prob);
	    double cost2 = ((h2[i] == -1) || (h2[i] == allele)) ? err_prob : (1.0 - err_prob);
	    if (i == 0) {
	      dp[i].score = std::min(std::log(cost1), std::log(cost2));
	      dp[i].hp = (std::log(cost1) < std::log(cost2)) ? 0 : 1;
	    } else {
	      double s1 = dp[i-1].score + std::log(cost1) + ((dp[i-1].hp == 0) ? 0 : switch_penalty);
	      double s2 = dp[i-1].score + std::log(cost2) + ((dp[i-1].hp == 1) ? 0 : switch_penalty);
	      dp[i].score = std::min(s1,s2);
	      dp[i].hp = (s1 < s2) ? 0 : 1;
	    }
          }
          new_partition[r_idx] = dp[block_size-1].hp;
          current_mec_score -= dp[block_size - 1].score;
          if (new_partition[r_idx] != partition[r_idx]) changed = true;
	}
        partition = new_partition;
        if (current_mec_score < min_mec_score) {
          min_mec_score = current_mec_score;
          best_partition = partition;
        }
        if (!changed) break;
      }
      
      // Build haplotypes
      std::vector<int> final_h1(block_size, -1);
      std::vector<int> final_h2(block_size, -1);
      for (int i = 0; i < block_size; ++i) {
	double h1_ref = 0;
	double h1_alt = 0;
	double h2_ref = 0;
	double h2_alt = 0;
	for(uint32_t r_idx = 0; r_idx < block_read_indices.size(); ++r_idx) {
	  ReadInfo const& ri = reads_map.at(read_hashes[block_read_indices[r_idx]]);
	  if (ri.mask[s+i]) {
	    double p_corr = 1.0 - _phredToErr(ri.bq[s+i]);
	    if (best_partition[r_idx] == 0) { // H1
	      if (ri.alt[s+i]) h1_alt += p_corr;
	      else h1_ref += p_corr;
	    } else { // H2
	      if (ri.alt[s+i]) h2_alt += p_corr;
	      else h2_ref += p_corr;
	    }
	  }
	}
	if (h1_ref + h1_alt > 0) final_h1[i] = (h1_alt > h1_ref) ? 1 : 0;
	if (h2_ref + h2_alt > 0) final_h2[i] = (h2_alt > h2_ref) ? 1 : 0;
      }
      for (int i = 0; i < block_size; ++i) {
	// Fill-in missing alleles
	if ((final_h1[i] != -1) && (final_h2[i] == -1)) final_h2[i] = 1 - final_h1[i];
	else if ((final_h1[i] == -1) && (final_h2[i] != -1)) final_h1[i] = 1 - final_h2[i];
	
	// Valid het. variant?
        if ((final_h1[i] != -1) && (final_h2[i] != -1) && (final_h1[i] != final_h2[i])) {
          out_hap[s + i] = final_h1[i];
          phase_block_id[s + i] = current_ps_id;
        }
      }
      ++current_ps_id;
    }
  }

  template<typename TConfig>
  inline void
  haplotagBam(TConfig const&c, std::vector<std::vector<Variant>> const& variants_by_chr) {
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    htsFile *out = sam_open(c.taggedfile.string().c_str(), "wb");
    if (sam_hdr_write(out, hdr) < 0) {
      std::cerr << "Warning: Could not write BAM header!" << std::endl;
    }

    // Variant positions
    std::vector<std::vector<int32_t>> posIndex(variants_by_chr.size());
    for (uint32_t refIndex = 0; refIndex < variants_by_chr.size(); ++refIndex) {
        const auto &vars = variants_by_chr[refIndex];
        posIndex[refIndex].reserve(vars.size());
        for (auto const &v : vars) posIndex[refIndex].push_back(v.pos);
    }

    uint32_t nTagged = 0;
    uint32_t nTotal = 0;
    bam1_t *rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      ++nTotal;
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) {
	if (sam_write1(out, hdr, rec) < 0) std::cerr << "Warning: Could not write BAM record!" << std::endl;
	continue;
      }
      if ((rec->core.tid < 0) || (rec->core.tid >= (int) variants_by_chr.size())) {
	if (sam_write1(out, hdr, rec) < 0) std::cerr << "Warning: Could not write BAM record!" << std::endl;
	continue;
      }

      const auto &vars = variants_by_chr[rec->core.tid];
      if (vars.empty()) {
	if (sam_write1(out, hdr, rec) < 0) std::cerr << "Warning: Could not write BAM record!" << std::endl;
	continue;
      }

      // Any variants for this read?
      auto lb = std::lower_bound(posIndex[rec->core.tid].begin(), posIndex[rec->core.tid].end(), rec->core.pos);
      auto ub = std::upper_bound(posIndex[rec->core.tid].begin(), posIndex[rec->core.tid].end(), bam_endpos(rec));
      if (lb == ub) {
	if (sam_write1(out, hdr, rec) < 0) std::cerr << "Warning: Could not write BAM record!" << std::endl;
	continue;
      }
      std::size_t sidx = lb - posIndex[rec->core.tid].begin();
      std::size_t eidx = ub - posIndex[rec->core.tid].begin();
      if (eidx > vars.size()) eidx = vars.size();

      // Load sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t const* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      
      // Count matches to each haplotype
      int hapCount[2] = {0, 0};
      int32_t dominantPS = 0;      
      for (std::size_t i = sidx; i < eidx; ++i) {
	const Variant &v = vars[i];
	if (v.hap == -1) continue;
	int read_offset = -1;

	// Parse cigar
	uint32_t const* cigar = bam_get_cigar(rec);
	int gp = rec->core.pos;
	int sp = 0;
	for (uint32_t k = 0; k < rec->core.n_cigar; ++k) {
	  int op = bam_cigar_op(cigar[k]);
	  int len = bam_cigar_oplen(cigar[k]);
	  if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
	    if ((v.pos >= gp) && (v.pos < gp + len)) {
	      read_offset = sp + (v.pos - gp);
	      break;
	    }
	    gp += len;
	    sp += len;
	  } else if (op == BAM_CINS) {
	    sp += len;
	  } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
	    gp += len;
	  } else if (op == BAM_CSOFT_CLIP) {
	    sp += len;
	  }
	}

	// Count HP support
	if ((read_offset >= 0) && (read_offset < rec->core.l_qseq)) {
	  std::string REF = v.ref;
	  std::string ALT = v.alt;
	  if (sequence.substr(read_offset, REF.size()) == REF) {
	    if (v.hap == 0) ++hapCount[0];
	    else ++hapCount[1];
	    if (dominantPS == 0) dominantPS = v.ps;
	  } else if (sequence.substr(read_offset, ALT.size()) == ALT) {
	    if (v.hap == 0) ++hapCount[1];
	    else ++hapCount[0];
	    if (dominantPS == 0) dominantPS = v.ps;
	  }
	}
      }

      // Tag read
      int hp = 0;
      if (hapCount[0] + hapCount[1] >= 2) {
	if (hapCount[0] >= 2 * hapCount[1]) hp = 1;
	else if (hapCount[1] >= 2 * hapCount[0]) hp = 2;
	if (hp) {
	  bam_aux_append(rec, "HP", 'i', 4, (uint8_t*)&hp);
	  bam_aux_append(rec, "PS", 'i', 4, (uint8_t*)&dominantPS);
	  ++nTagged;
	}
      }
      if (sam_write1(out, hdr, rec) < 0) std::cerr << "Warning: Could not write BAM record!" << std::endl;
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);
    sam_close(out);

    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Haplotagged " << nTagged << "/" << nTotal << " reads." << std::endl;
}


  template<typename TConfig>
  inline void
  writeVcf(TConfig const& c, const std::vector<std::vector<Variant>>& variants_by_chr) {

    // Load BAM file
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load BCF file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);
    bcf_hdr_append(bcfhdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">");
    
    // Fetch sample
    int32_t sampleIndex = -1;
    for (int32_t i = 0; i < bcf_hdr_nsamples(bcfhdr); ++i) {
      if (bcfhdr->samples[i] == c.sampleName) sampleIndex = i;
    }
    if (sampleIndex < 0) {
      if (bcf_hdr_nsamples(bcfhdr) > 0) {
	std::cerr << "Warning: None of the samples in the VCF/BCF file matches " << c.sampleName << "! Taking first sample in VCF/BCF for phasing: " << bcfhdr->samples[0] << std::endl;
	sampleIndex = 0;
      } else {
	std::cerr << "Warning: No genotypes are present in the VCF/BCF file. All variants are considered as het. bi-allelic variants." << std::endl;
      }	
    }

    // Output file
    std::string fmtout = "wb";
    if (c.outfile.string() == "-") fmtout = "w";
    htsFile* out_vcf = hts_open(c.outfile.string().c_str(), fmtout.c_str());
    if (bcf_hdr_write(out_vcf, bcfhdr) < 0) std::cerr << "Warning: Could not write BCF header!" << std::endl;

    // Iterate input VCF
    int32_t lastRefIdx = -1;
    int32_t refIndex = -1;
    uint32_t ptr = 0;
    bcf1_t* rec = bcf_init();
    while (bcf_read(ibcffile, bcfhdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_ALL);

      // New chromosome?
      if (rec->rid != lastRefIdx) {
	// Match chr rid to BAM tid
	std::string chrom = bcf_hdr_id2name(bcfhdr, rec->rid);
	refIndex = bam_name2id(hdr, chrom.c_str());
	ptr = 0;
	lastRefIdx = rec->rid;
      }
      if (refIndex < 0) {
	if (bcf_write(out_vcf, bcfhdr, rec) < 0) std::cerr << "Warning: Could not write BCF record!" << std::endl;
	continue;
      }

      // Find variant
      while ((ptr < variants_by_chr[refIndex].size()) && (variants_by_chr[refIndex][ptr].pos < rec->pos)) ++ptr;

      // Update PS and GT
      if ((ptr < variants_by_chr[refIndex].size()) && (variants_by_chr[refIndex][ptr].pos == rec->pos) && (rec->n_allele == 2) && (std::string(rec->d.allele[0]) == variants_by_chr[refIndex][ptr].ref) && (std::string(rec->d.allele[1]) == variants_by_chr[refIndex][ptr].alt)) {
	int32_t ngt_arr=0;
	int32_t *gt_arr=NULL;
	bcf_get_genotypes(bcfhdr, rec, &gt_arr, &ngt_arr);
	int hap = variants_by_chr[refIndex][ptr].hap;
	if (hap == -1) {
	  gt_arr[sampleIndex*2] = bcf_gt_unphased(0);
	  gt_arr[sampleIndex*2 + 1] = bcf_gt_unphased(1);
	  bcf_update_format_int32(bcfhdr, rec, "PS", NULL, 0);
	} else {
	  gt_arr[sampleIndex*2] = bcf_gt_phased(hap ? 1 : 0);
	  gt_arr[sampleIndex*2 + 1] = bcf_gt_phased(hap ? 0 : 1);
	  int32_t psVal = variants_by_chr[refIndex][ptr].ps;
	  bcf_update_format_int32(bcfhdr, rec, "PS", &psVal, 1);	  
	}
	bcf_update_genotypes(bcfhdr, rec, gt_arr, ngt_arr);
	if (gt_arr) free(gt_arr);
      }
      if (bcf_write(out_vcf, bcfhdr, rec) < 0) std::cerr << "Warning: Could not write BCF record!" << std::endl;
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(bcfhdr);
    bcf_close(ibcffile);
    // Close output file
    bcf_close(out_vcf);
    if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);
  }

  template<typename TConfig>
  inline int
  runPhase(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("cameo.prof");
#endif
    
    // Parse VCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse VCF/BCF variants" << std::endl;
    std::vector<std::vector<Variant>> variants_by_chr(c.nchr);
    parseVcfVariants(c, variants_by_chr);

    // Iterate chromosomes
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Phase variants" << std::endl;
    ThreadPool pool(c.maxThreads);
    for (int refIndex = 0; refIndex < c.nchr; ++refIndex) {
      auto &vars = variants_by_chr[refIndex];
      if (vars.empty()) continue;

      // Global genomic position to variant index map for this chr
      std::unordered_map<int,int> global_pos_to_idx;
      global_pos_to_idx.reserve(vars.size()*2);
      for (int vi = 0; vi < (int)vars.size(); ++vi) global_pos_to_idx[vars[vi].pos] = vi;
      
      // Store results per chunk
      int chrom_len_est = vars.back().pos + 1;
      int nchunks = (chrom_len_est + c.chunkSize - 1) / c.chunkSize;	
      std::vector<std::shared_ptr<std::vector<int>>> chunk_haps(nchunks);
      std::vector<std::shared_ptr<std::vector<int>>> chunk_positions(nchunks);
      std::vector<std::shared_ptr<std::vector<int>>> chunk_psids(nchunks);
      
      // Enqueue chunk tasks
      std::vector<std::future<void>> futures;
      for (int chunk_id = 0; chunk_id < nchunks; ++chunk_id) {
	int chunk_start = chunk_id * c.chunkSize;
	int chunk_end = std::min(chunk_start + c.chunkSize + c.chunkOverlap, chrom_len_est);
	
	// Build chunk_vars and pos_to_local
	std::vector<Variant> chunk_vars;
	chunk_vars.reserve(4096);
	std::unordered_map<int,int> pos_to_local;
	pos_to_local.reserve(4096);
	for (int vi = 0; vi < (int)vars.size(); ++vi) {
	  int p = vars[vi].pos;
	  if ((p >= chunk_start) && (p < chunk_end)) {
	    pos_to_local[p] = (int)chunk_vars.size();
	    chunk_vars.push_back(vars[vi]);
	  }
	}
	if (chunk_vars.empty()) continue;
	
	// Collect reads for each chunk and phase
	chunk_haps[chunk_id] = std::make_shared<std::vector<int>>();
	chunk_positions[chunk_id] = std::make_shared<std::vector<int>>();
	chunk_psids[chunk_id] = std::make_shared<std::vector<int>>();
	futures.push_back(pool.enqueue([=, &chunk_haps, &chunk_positions, &chunk_psids]() {
	  std::unordered_map<std::size_t, ReadInfo> reads_map;
	  collectReadsPerChunk(c, refIndex, chunk_vars, pos_to_local, reads_map, chunk_start, chunk_end);
	  
	  // Solve chunk
	  std::vector<int> hap;
	  std::vector<int> psids;
	  solveChunk(c, chunk_vars, reads_map, hap, psids);
	  
	  auto hp = chunk_haps[chunk_id];
	  auto posv = chunk_positions[chunk_id];
	  auto psv = chunk_psids[chunk_id];
	  hp->assign(hap.begin(), hap.end());
	  posv->reserve(chunk_vars.size());
	  psv->assign(psids.begin(), psids.end());
	  for (uint32_t i = 0; i < chunk_vars.size(); ++i) posv->push_back(chunk_vars[i].pos);
	}));
      }
      pool.waitAll();

      // Merge
      for (int chunk_id = 0; chunk_id < nchunks; ++chunk_id) {
        auto const& hp = chunk_haps[chunk_id];
        auto const& posv = chunk_positions[chunk_id];
        auto const& psv = chunk_psids[chunk_id];
        if ((!hp) || (hp->empty())) continue;
        int chunk_start = chunk_id * c.chunkSize;
        
        // Stitch on chunk overlap?
        std::map<int32_t, int> ps_vote;
        std::map<int32_t, int32_t> ps_target;
        if (chunk_id > 0) {
          int overlap_start = chunk_start;
          int overlap_end = chunk_start + c.chunkOverlap;
          for(uint32_t k = 0; k < posv->size(); ++k) {
            int gpos = (*posv)[k];
            if (gpos >= overlap_start && gpos < overlap_end) {
              auto it = global_pos_to_idx.find(gpos);
              if (it != global_pos_to_idx.end()) {
                int gidx = it->second;
                if (vars[gidx].hap != -1 && (*hp)[k] != -1 && vars[gidx].ps != 0) {
                  ps_target[(*psv)[k]] = vars[gidx].ps;
                  if (vars[gidx].hap == (*hp)[k]) ++ps_vote[(*psv)[k]];
                  else --ps_vote[(*psv)[k]];
                }
              }
            }
          }
        }
        for(uint32_t k = 0; k < posv->size(); ++k) {
          int gpos = (*posv)[k];
          //if ((chunk_id > 0) && (gpos < chunk_start)) continue;
          auto it = global_pos_to_idx.find(gpos);
          if (it != global_pos_to_idx.end()) {
	    if(vars[it->second].ps != 0) continue;  // Skip already phased variants from the previous chunk
            int gidx = it->second;
            int32_t current_ps = (*psv)[k];
            int current_hap = (*hp)[k];

            if ((current_hap != -1) && (current_ps != 0)) {
              auto vote_it = ps_vote.find(current_ps);
              if ((vote_it != ps_vote.end()) && (std::abs(vote_it->second) >= 1)) {
		// Stitch
                vars[gidx].ps = ps_target[current_ps];
                vars[gidx].hap = (vote_it->second > 0) ? current_hap : 1 - current_hap;
              } else {
		// No stitching
                vars[gidx].ps = (chunk_id << 20) | current_ps; // Create a unique tmp PS by shifting 20 bits
                vars[gidx].hap = current_hap;
              }
            } else {
              vars[gidx].hap = -1;
              vars[gidx].ps = 0;
            }
          }
        }
      }
      // Remap PS Identifiers
      std::map<int32_t, int32_t> final_ps_map;
      int32_t global_ps_counter = 1;
      for (auto& var : vars) {
        if (var.ps != 0) {
          if (final_ps_map.find(var.ps) == final_ps_map.end()) final_ps_map[var.ps] = global_ps_counter++;
          var.ps = final_ps_map[var.ps];
        }
      }

      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Chromosome index " << (refIndex + 1) << " done (variants=" << vars.size() << ", chunks=" << nchunks << ")" << std::endl;
    }


    // Write VCF
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Write phased BCF" << std::endl;
    writeVcf(c, variants_by_chr);

    // Haplotag BAM
    if (c.hasTaggedBam) {
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Write tagged BAM" << std::endl;
      haplotagBam(c, variants_by_chr);
    }

#ifdef PROFILE
    ProfilerStop();
#endif

    // Done
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    
    return 0;
  }


  int phase(int argc, char **argv) {
    PhaseConfig c;
   
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("vcffile,f", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "VCF/BCF file with variants")
      ("base-qual,b", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(1), "min. base quality")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName), "sample name")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output phased BCF file")
      ("tagged,d", boost::program_options::value<boost::filesystem::path>(&c.taggedfile), "output tagged BAM file [optional]")
      ("threads,t", boost::program_options::value<uint32_t>(&c.maxThreads)->default_value(24), "number of threads")
     ;
    
    boost::program_options::options_description phasing("Phasing options");
    phasing.add_options()
      ("chunksize,u", boost::program_options::value<int32_t>(&c.chunkSize)->default_value(1000000), "chunk size")
      ("chunkoverlap,l", boost::program_options::value<int32_t>(&c.chunkOverlap)->default_value(50000), "chunk overlap")
      ("maxiter,m", boost::program_options::value<uint32_t>(&c.maxIter)->default_value(50), "max. iterations per block")
      ("minreadsperblock,n", boost::program_options::value<uint32_t>(&c.minReadsPerBlock)->default_value(2), "min. reads per block")
      ("maxreadsperchunk,r", boost::program_options::value<uint32_t>(&c.maxReadsPerChunk)->default_value(10000), "max. reads per chunk")
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
      std::cerr << "Usage: cameo " << argv[0] << " [OPTIONS] -g <ref.fa> -f <input.vcf> <input.bam>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }

    // Tagged BAM?
    if (vm.count("tagged")) c.hasTaggedBam = true;
    else c.hasTaggedBam = false;
    
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
   
    return runPhase(c);
  }
  
}

#endif
