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

namespace cameo {
  
  struct PhaseConfig {
    bool hasTaggedBam;
    uint16_t minBaseQual;
    uint16_t minMapQual;
    uint32_t maxThreads;
    uint32_t minOverlapVar;
    uint32_t maxReadsPerChunk;
    uint32_t maxIter;
    int32_t nchr;
    int32_t minCov;
    int32_t chunkSize;
    int32_t chunkOverlap;
    int32_t minReadsPerBlock;
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
    explicit Variant(int32_t const p) : pos(p), ps(0), hap(0), ref(""), alt("") {}
    Variant(int32_t const p, std::string const& r, std::string const& a) : pos(p), ps(0), hap(0), ref(r), alt(a) {}

    bool operator<(const Variant& v2) const {
      return (pos<v2.pos);
    }
  };

  struct ReadInfo {
    boost::dynamic_bitset<> mask;
    boost::dynamic_bitset<> alt;
    std::vector<uint8_t> bq;
    uint16_t mapq;
    
    ReadInfo() : mask(), alt(), bq() {}
    explicit ReadInfo(int const numVars) {
      mask.resize(numVars, false);
      alt.resize(numVars, false);
      bq.resize(numVars, 0);
    }    
  };

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
	  if (rAllele != boost::to_upper_copy(std::string(seq + rec->pos, seq + rec->pos + rAllele.size()))) {
	    std::cerr << "Warning: REF allele mismatch at " << hdr->target_name[refIndex] << ":" << (rec->pos + 1) << std::endl;
	  } else {
	    variants_by_chr[refIndex].push_back(Variant(rec->pos, std::string(rec->d.allele[0]), std::string(rec->d.allele[1])));
	    ++kept;
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
    for(uint32_t refIndex = 0; refIndex < variants_by_chr.size(); ++refIndex) {
      std::sort(variants_by_chr[refIndex].begin(), variants_by_chr[refIndex].end());
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


      hts_itr_t* itr = sam_itr_queryi(idx, refIndex, chunk_genomic_start, chunk_genomic_end);
      bam1_t* rec = bam_init1();
      int numVars = chunk_vars.size();
      while (sam_itr_next(samfile, itr, rec) >= 0) {
        if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
        if ((rec->core.l_qseq <= 0) || (rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	// Get base qualities
	typedef std::vector<uint8_t> TQuality;
	TQuality quality;
	quality.resize(rec->core.l_qseq);
	uint8_t* qualptr = bam_get_qual(rec);
	for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
	
	// Add ReadInfo
	std::size_t seed = hash_lr(bam_get_qname(rec));
	readTmp[seed] = ReadInfo(numVars);
	ReadInfo &ri = readTmp[seed];

	// Parse cigar
        int gp = rec->core.pos;
        int sp = 0;
	uint32_t* cigar = bam_get_cigar(rec);
        for (uint32_t ci = 0; ci < rec->core.n_cigar; ++ci) {
	  if ((bam_cigar_op(cigar[ci]) == BAM_CMATCH) || (bam_cigar_op(cigar[ci]) == BAM_CEQUAL) || (bam_cigar_op(cigar[ci]) == BAM_CDIFF)) {
	    for (uint32_t i = 0; i < bam_cigar_oplen(cigar[ci]); ++i) {
	      int rp = gp + i;
	      if ((rp < chunk_genomic_start) || (rp >= chunk_genomic_end)) continue;
	      auto itpos = pos_to_local.find(rp);
	      if (itpos == pos_to_local.end()) continue;
	      int local = itpos->second;
	      int qi = sp + i;
	      if ((qi < 0) || (qi >= rec->core.l_qseq)) continue;
	      if (quality[qi] < c.minBaseQual) continue;
	      std::string REF = chunk_vars[local].ref;
	      std::string ALT = chunk_vars[local].alt;
	      if (sequence.substr(qi, REF.size()) == REF) {
		ri.mask[local] = true;
		ri.bq[local] = quality[qi];
		ri.mapq = rec->core.qual;
	      }
	      else if (sequence.substr(qi, ALT.size()) == ALT) {
		ri.mask[local] = true;
		ri.alt[local] = true;
		ri.bq[local] = quality[qi];
		ri.mapq = rec->core.qual;
	      }
	    }
	    gp += bam_cigar_oplen(cigar[ci]);
	    sp += bam_cigar_oplen(cigar[ci]);
	  } else if (bam_cigar_op(cigar[ci]) == BAM_CINS) {
	    sp += bam_cigar_oplen(cigar[ci]);
	  } else if ((bam_cigar_op(cigar[ci]) == BAM_CDEL) || (bam_cigar_op(cigar[ci]) == BAM_CREF_SKIP)) {
	    gp += bam_cigar_oplen(cigar[ci]);
	  } else if (bam_cigar_op(cigar[ci]) == BAM_CSOFT_CLIP) {
	    sp += bam_cigar_oplen(cigar[ci]);
	  }
        }
      }
      bam_destroy1(rec);
      hts_itr_destroy(itr);
      hts_idx_destroy(idx);
      bam_hdr_destroy(hdr);
      sam_close(samfile);

      // Remove reads that do not cover enough variants and downsample (if needed)
      std::vector<std::pair<uint32_t, std::size_t> > varCounts;
      for(typename TReadMap::iterator it = readTmp.begin(); it != readTmp.end(); ++it) varCounts.push_back(std::make_pair(it->second.mask.count(), it->first));
      std::sort(varCounts.begin(), varCounts.end(), std::greater<>());
      uint32_t rcount = 0;
      for(uint32_t i = 0; i < varCounts.size(); ++i) {
	if (varCounts[i].first >= c.minOverlapVar) {
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
    int B = chunk_vars.size();
    out_hap.assign(B, 0);
    phase_block_id.assign(B, 0);

    int R = reads_map.size();
    if (R < c.minReadsPerBlock) return;

    // Collect read masks and alt bitsets
    std::vector<boost::dynamic_bitset<>> read_mask;
    std::vector<boost::dynamic_bitset<>> read_alt;
    read_mask.reserve(R);
    read_alt.reserve(R);
    for (auto const &kv : reads_map) {
      read_mask.push_back(kv.second.mask);
      read_alt.push_back(kv.second.alt);
    }

    // Read span count between consecutive variants
    std::vector<int> spanCount(std::max(0, B - 1), 0);
    for (int r = 0; r < R; ++r) {
      int last = -1;
      for (int i = 0; i < B; ++i) {
	if (read_mask[r][i]) {
	  if (last >= 0) ++spanCount[last];
	  last = i;
	}
      }
    }

    // Split into phased blocks
    std::vector<std::pair<int,int>> blocks;
    int block_start = 0;
    for (int i = 0; i < B - 1; ++i) {
      if (spanCount[i] < c.minReadsPerBlock) {
	blocks.push_back(std::make_pair(block_start, i));
	block_start = i + 1;
      }
    }
    blocks.push_back(std::make_pair(block_start, B - 1));

    // Process blocks
    int global_ps_id = 1;
    for (auto const &blk : blocks) {
      int s = blk.first;
      int e = blk.second;
      int block_size = e - s + 1;
      if (block_size <= 0) continue;
      
      // If only one variant in this block
      if (block_size == 1) {
	out_hap[s] = 0;
	phase_block_id[s] = global_ps_id++;
	continue;
      }
      
      // Weighted adjacency between variant i and j
      std::vector<std::unordered_map<int,double>> adj(block_size);
      for (int r = 0; r < R; ++r) {
	std::vector<int> idxs;
	std::vector<int> alts;
	idxs.reserve(block_size);
	alts.reserve(block_size);
	for (int local = 0; local < block_size; ++local) {
	  int gidx = s + local;
	  if (read_mask[r][gidx]) {
	    idxs.push_back(local);
	    alts.push_back(read_alt[r][gidx] ? 1 : 0);
	  }
	}
	int k = idxs.size();
	if (k < 2) continue;
	
	// Calculate edge weight
	double span = std::sqrt((double)k);
	double w = std::min(4.0, std::max(0.5, span));
	for (int i = 0; i < k; ++i) {
	  int vi = idxs[i];
	  int ai = alts[i];
	  for (int j = i + 1; j < k; ++j) {
	    int vj = idxs[j];
	    int aj = alts[j];
	    double delta = (ai == aj) ? w : -w;
	    adj[vi][vj] += delta;
	    adj[vj][vi] += delta;
	  }
	}
      }
      
      // No spanning reads?
      bool hasEdges = false;
      for (int i = 0; i < block_size; ++i) {
	if (!adj[i].empty()) {
	  hasEdges = true;
	  break;
	}
      }
      if (!hasEdges) {
	for (int local = 0; local < block_size; ++local) {
	  int cov = 0, altc = 0;
	  int gidx = s + local;
	  for (int r = 0; r < R; ++r) {
	    if (read_mask[r][gidx]) {
	      ++cov;
	      if (read_alt[r][gidx]) ++altc;
	    }
	  }
	  out_hap[gidx] = (cov > 0 && altc * 2 >= cov) ? 1 : 0;
	  phase_block_id[gidx] = global_ps_id;
	}
	++global_ps_id;
	continue;
      }

      // Initialize haplotype signs using per-site majority
      std::vector<int> hap_sign(block_size, 1);
      for (int local = 0; local < block_size; ++local) {
	int cov = 0;
	int altc = 0;
	int gidx = s + local;
	for (int r = 0; r < R; ++r) {
	  if (read_mask[r][gidx]) {
	    ++cov;
	    if (read_alt[r][gidx]) ++altc;
	  }
	}
	hap_sign[local] = ((cov > 0) && (altc * 2 >= cov)) ? 1 : -1;
      }
      
      // Iterative linkage propagation
      const int MAX_ITERS = std::max(5, (int)c.maxIter);
      bool converged = false;
      for (int iter = 0; iter < MAX_ITERS && !converged; ++iter) {
	converged = true;
	std::vector<int> new_sign(block_size, 1);
	for (int i = 0; i < block_size; ++i) {
	  double score = 0.0;
	  for (auto const &kv : adj[i]) {
	    score += kv.second * hap_sign[kv.first];
	  }
	  int ns = (score > 0.0) ? 1 : ((score < 0.0) ? -1 : hap_sign[i]);
	  if (ns != hap_sign[i]) converged = false;
	  new_sign[i] = ns;
	}
	hap_sign.swap(new_sign);
      }

      // Read reassignment to haplotypes
      std::vector<int> read_assign(R, -1);
      for (int r = 0; r < R; ++r) {
	double score = 0.0;
	for (int local = 0; local < block_size; ++local) {
	  int gidx = s + local;
	  if (!read_mask[r][gidx]) continue;
	  int readAllele = read_alt[r][gidx] ? 1 : 0;
	  int hapAllele = (hap_sign[local] == 1) ? 1 : 0;
	  score += (readAllele == hapAllele) ? 1.0 : -1.0;
	}
	read_assign[r] = (score >= 0.0) ? 0 : 1;
      }

      // Weighting
      for (int local = 0; local < block_size; ++local) {
	double vote = 0.0, weight_sum = 0.0;
	int gidx = s + local;
	for (int r = 0; r < R; ++r) {
	  if (!read_mask[r][gidx]) continue;
	  int readAllele = read_alt[r][gidx] ? 1 : 0;
	  int hapAllele = (hap_sign[local] == 1) ? 1 : 0;
	  int contribution = (readAllele == hapAllele) ? 1 : -1;
	  if (read_assign[r] == 1) contribution = -contribution;
	  double w = 1.0; // MAPQ/BQ weights?
	  vote += w * contribution;
	  weight_sum += w;
	}
	
	// Compute phasing confidence
	double conf = (weight_sum > 0.0) ? std::abs(vote) / weight_sum : 0.0;
	double linkage_strength = 0.0;
	for (auto const& kv : adj[local]) linkage_strength += std::abs(kv.second);
	linkage_strength /= std::max(1.0, (double)adj[local].size());
	double combined_conf = 0.5 * conf + 0.5 * std::min(1.0, linkage_strength / 5.0);
	
	// Output haplotype
	double thr = 0.45 + 0.15 * std::tanh((weight_sum - 8.0) / 6.0);
	if (combined_conf < thr) out_hap[gidx] = -1;  // unphased
	else {
	  out_hap[gidx] = (vote >= 0.0) ? ((hap_sign[local] == 1) ? 1 : 0) : ((hap_sign[local] == 1) ? 0 : 1);
	}
	phase_block_id[gidx] = global_ps_id;
      }
      ++global_ps_id;
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
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      
      // Count matches to each haplotype
      int hapCount[2] = {0, 0};
      int32_t dominantPS = 0;      
      for (std::size_t i = sidx; i < eidx; ++i) {
	const Variant &v = vars[i];
	int read_offset = -1;

	// Parse cigar
	uint32_t *cigar = bam_get_cigar(rec);
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
	    dominantPS = vars[i].ps;
	  } else if (sequence.substr(read_offset, ALT.size()) == ALT) {
	    if (v.hap == 0) ++hapCount[1];
	    else ++hapCount[0];
	    dominantPS = vars[i].ps;
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
    int32_t lastpos = -1;
    int32_t refIndex = -1;
    uint32_t ptr = 0;
    bcf1_t* rec = bcf_init();
    while (bcf_read(ibcffile, bcfhdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_STR);

      // New chromosome?
      if (rec->rid != lastRefIdx) {
	// Match chr rid to BAM tid
	std::string chrom = bcf_hdr_id2name(bcfhdr, rec->rid);
	refIndex = bam_name2id(hdr, chrom.c_str());
	if (refIndex < 0) {
	  std::cerr << "Warning: Chromosome does not exist in BAM file: " << chrom << std::endl;
	  continue;
	}
	ptr = 0;
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
	if ((a0 == a1) || (a0 < 0) || (a1 < 0) || (!((a0==0 && a1==1) || (a0==1 && a1==0)))) {
	  if (gt_arr) free(gt_arr);
	  continue;
	}
	if (rec->pos != lastpos) {
	  if ((ptr < variants_by_chr[refIndex].size()) && (variants_by_chr[refIndex][ptr].pos == rec->pos)) {
	    int hap = variants_by_chr[refIndex][ptr].hap;
	    if (hap == -1) {
	      gt_arr[sampleIndex*2] = bcf_gt_unphased(0);
	      gt_arr[sampleIndex*2 + 1] = bcf_gt_unphased(1);
	    } else {
	      gt_arr[sampleIndex*2] = bcf_gt_phased(hap ? 1 : 0);
	      gt_arr[sampleIndex*2 + 1] = bcf_gt_phased(hap ? 0 : 1);
	    }
	    bcf_update_genotypes(bcfhdr, rec, gt_arr, ngt_arr);
	    int32_t psVal = variants_by_chr[refIndex][ptr].ps;
	    if (hap != -1) bcf_update_format_int32(bcfhdr, rec, "PS", &psVal, 1);
	    else bcf_update_format_int32(bcfhdr, rec, "PS", NULL, 0);
	    ++ptr;
	  } else {
	    std::cerr << "Warning: Variants are out of sync between parsing and writing final VCF!" << std::endl;
	  }
	  if (bcf_write(out_vcf, bcfhdr, rec) < 0) std::cerr << "Warning: Could not write BCF record!" << std::endl;
	  lastpos = rec->pos;
	}
	if (gt_arr) free(gt_arr);
      }
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
    uint32_t lastChrPsVal = 0;
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
	  for (size_t i = 0; i < chunk_vars.size(); ++i) posv->push_back(chunk_vars[i].pos);
	}));
      }
      pool.waitAll();

      // Merge chunk
      for (int chunk_id = 0; chunk_id < nchunks; ++chunk_id) {
	auto hp = chunk_haps[chunk_id];
	auto posv = chunk_positions[chunk_id];
	auto psv = chunk_psids[chunk_id];
	if (!hp || !posv || !psv) continue;
	for (size_t k = 0; k < posv->size(); ++k) {
	  int gpos = (*posv)[k];
	  auto it = global_pos_to_idx.find(gpos);
	  if (it != global_pos_to_idx.end()) {
	    int gidx = it->second;
	    variants_by_chr[refIndex][gidx].hap = (*hp)[k];
	    variants_by_chr[refIndex][gidx].ps = (*psv)[k];
	  }
	}
	hp->clear();
	posv->clear();
	psv->clear();
      }

      // Update genome-wide PS values
      if (!variants_by_chr[refIndex].empty()) {
	uint32_t newChrPsVal = lastChrPsVal;
	for(uint32_t i = 0; i < variants_by_chr[refIndex].size(); ++i) {
	  if (variants_by_chr[refIndex][i].hap != -1) {
	    variants_by_chr[refIndex][i].ps += lastChrPsVal;
	    newChrPsVal = variants_by_chr[refIndex][i].ps;
	  }
	}
	lastChrPsVal = newChrPsVal;
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
      ("base-qual,b", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(10), "min. base quality")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName), "sample name")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output phased BCF file")
      ("tagged,d", boost::program_options::value<boost::filesystem::path>(&c.taggedfile), "output tagged BAM file [optional]")
      ("threads,t", boost::program_options::value<uint32_t>(&c.maxThreads)->default_value(24), "number of threads")
     ;
    
    boost::program_options::options_description phasing("Phasing options");
    phasing.add_options()
      ("mincov,c", boost::program_options::value<int32_t>(&c.minCov)->default_value(10), "min. variant coverage")
      ("chunksize,u", boost::program_options::value<int32_t>(&c.chunkSize)->default_value(1000000), "chunk size")
      ("chunkoverlap,l", boost::program_options::value<int32_t>(&c.chunkOverlap)->default_value(50000), "chunk overlap")
      ("minoverlap,p", boost::program_options::value<uint32_t>(&c.minOverlapVar)->default_value(3), "min. variants spanned by a read")
      ("maxiter,m", boost::program_options::value<uint32_t>(&c.maxIter)->default_value(50), "max. iterations per block")
      ("minreadsperblock,n", boost::program_options::value<int32_t>(&c.minReadsPerBlock)->default_value(2), "min. reads per block")
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
