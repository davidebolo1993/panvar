#pragma once

// Recruiting reads to locus indexes without aligning: each read (or pair) is screened by its closed
// syncmers against the loaded loci; for every locus it hits, each of its k-mers is labelled P, S, R (the
// locus table), D (the decoys) or N (neither), each read's locus segment is found from the labels, and the
// evidence in the segments decides the fragment's fate: certain, ambiguous (with a reason) or unrecruited.

#include "panvar/locus_index.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace panvar::recruit {

enum class Method { Runs, Kmers, Fraction };

struct Evidence {
    // Within each read's locus segment (the stretch of its k-mers that reads as locus sequence):
    // independent sites, private (or decoy) k-mers counted only when they start at least k bases after the
    // last one counted, so all k-mers a single differing base touches count once; a k-mer already counted
    // in the other mate is not counted again.
    std::uint32_t sites_p = 0;
    std::uint32_t sites_d = 0;
    std::uint32_t decoy_out = 0;   // decoy sites outside the segments
    std::uint32_t runs_n = 0;      // runs of k-mers in neither table (a sequencing error or a new allele each)
    std::uint32_t seg_kmers = 0;   // callable k-mers in the segments
    // Whole reads:
    std::uint32_t n_p = 0, n_s = 0, n_r = 0, n_d = 0, n_n = 0;   // k-mers per label
    std::uint32_t kmers = 0;       // k-mers read (windows without N bases)
    std::uint32_t windows = 0;     // k-mer windows the read has, callable or not
    std::uint32_t bases = 0;
    std::uint32_t multi = 0;       // P k-mers that are multi-copy within a haplotype
    void add(const Evidence& o);
};

struct ReadOptions {
    std::vector<std::string> indexes;    // .plx files
    std::string reads1, reads2;          // FASTQ (plain or gzip); reads2 for pairs
    std::string out;                     // per-fragment table (TSV)
    Method method = Method::Runs;
    std::uint32_t min_sites = 2;         // private sites for certain
    double decoy_share = 0.05;           // decoy sites tolerated, as a share of the private sites (runs method)
    double min_density = 0.05;           // private sites per 31 callable bases: chance hits in a long read are sparse
    double error_rate = 0.002;           // per base: the expected mismatch runs, beside one new allele per kb
    double min_fraction = 0.5;           // the fraction method's cut
    std::size_t threads = 2;
    std::function<void(const std::string&)> log;
};

struct ReadStats {
    std::uint64_t fragments = 0, screened = 0;
    std::uint64_t certain = 0, ambiguous = 0, unrecruited = 0;
};

ReadStats recruit_reads(const ReadOptions& options);

}  // namespace panvar::recruit
