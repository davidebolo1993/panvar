#pragma once

// The locus index: every canonical 31-mer of a locus graph's haplotypes, labelled against the background.
//   P (private)  outside the locus, the background holds it nowhere
//   S (shared)   1 to --informative-max times elsewhere per background haplotype (a FASTA genome is one,
//                a graph genome one per sample#haplotype of its paths), on average
//   R (repeat)   more often than that
// "Outside" is the background count minus the locus's own copies in every background genome, found
// through the anchors and read from each genome's FASTA (the coordinate reference: its exact interval).
// Beside the labels: the graph nodes each k-mer starts on, its largest copy number within one
// haplotype, the backbone flag (single-copy, in every haplotype), the decoy k-mers (the homologous
// windows' own alleles), the screen (closed syncmers among P and S k-mers) and per-length ceilings.
// Evidence is counted in independent sites: a run of r consecutive private k-mers holds ceil(r / k) of
// them, since one differing base changes at most k consecutive k-mers.

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace panvar::recruit {

constexpr std::uint8_t kPrivate = 0, kShared = 1, kRepeat = 2;

struct GenomeSpan {
    std::string genome, contig;
    std::uint64_t start = 0, end = 0;
    bool forward = true;
    std::uint64_t subtracted = 0;      // locus k-mer occurrences found in it
    bool piece = false;                // a piece of the locus off the located copies (no anchors placed it)
    std::uint64_t anchored_start = 0, anchored_end = 0;   // between the outermost inside anchors
    std::uint64_t locus_start = 0, locus_end = 0;         // that stretch walked out through the locus's k-mers: the
                                                          // copy; the rest of [start, end) is its uncertain ends
};

struct DecoyWindow {
    std::string genome, contig;
    std::uint64_t start = 0, end = 0;
    std::uint64_t hits = 0;            // shared k-mer hits that made it
};

struct Ceiling {
    std::size_t length = 0;            // read length
    double p3_min = 0, p3_median = 0;  // fraction of windows with >= 3 private k-mers (min, median over paths)
    double sites2_min = 0, sites2_median = 0;        // with >= 2 independent private sites
    double single2_min = 0, single2_median = 0;      // with >= 2 independent private single-copy sites
};

// Calibration of one read class: reads simulated from the locus haplotypes of held-out samples (their
// own k-mers removed from the index), from the decoy windows and from the coordinate reference's flanks,
// scored as recruit scores them; the counts certain at each private-site cut m (index m, 1 to kMaxCut).
constexpr std::uint32_t kMaxCut = 8;
// Scores are binned at this width; a fragment reaches a tier when its score's bin is at least the tier's.
constexpr double kScoreBin = 1.0 / 16;
struct ScoreHistogram {
    std::vector<std::int32_t> bin;            // increasing
    std::vector<std::uint64_t> count;
    std::uint64_t at_least(std::int32_t b) const;   // fragments in bins >= b
};
// A calibration genome's off-target fragments in one read class.
struct GenomeTally {
    std::string genome;
    std::uint64_t fragments = 0, screened = 0;
    std::uint64_t with_private = 0;           // screened with a private site under its relabelling
    std::uint64_t copies = 0;                 // fragments left out: over a stretch that reads as a copy of the locus
    double r_out = 0, q_out = 0;              // its own off-target site rates
};
struct Calibration {
    std::string read_class;                   // short, hifi, ont
    double error_rate = 0;
    std::uint64_t read_length = 0, insert = 0;   // insert 0: single reads
    std::uint64_t held_out = 0;               // samples held out
    std::uint64_t locus = 0;                  // their fragments
    std::uint64_t flank = 0, decoy = 0;       // off-target fragments simulated
    std::uint64_t flank_screened = 0, decoy_screened = 0;
    std::array<std::uint64_t, kMaxCut + 1> locus_certain{}, flank_certain{}, decoy_certain{};
    // The score: private and decoy sites per k bases of judged sequence in held-out locus fragments (r_in,
    // q_in) and in off-target ones (r_out, q_out: the median over the calibration genomes of each one's own rates,
    // so no single genome sets them; without calibration genomes, the simulated off-target fragments'), all
    // screened; and histograms of the deciding score (bins of kScoreBin) of the fragments that can be called,
    // held-out, off-target simulated from the flanks and decoy windows, and off-target scanned from the
    // calibration genomes (each held out of the background in turn).
    double r_in = 0, r_out = 0, q_in = 0, q_out = 0;
    std::uint64_t off_screened = 0;           // the simulated off-target fragments the score screened (its long reads
                                              // over the decoy windows are laid apart from the first rules')
    std::uint64_t genome_fragments = 0, genome_screened = 0;
    ScoreHistogram locus_score, off_score, genome_score;
    std::vector<GenomeTally> genomes;         // per calibration genome
};

// A genome held out of the background in calibration: its homologous windows (with context) and the flanks of its
// copies of the locus are read as a new sample's off-target sequence, with the locus k-mers only it holds outside
// the locus labelled private and the decoy k-mers only it holds unknown.
struct CalibrationGenome {
    std::string genome;
    std::uint64_t windows = 0, bases = 0;      // homologous windows, and bases read (windows, context, flanks)
    std::uint64_t homolog_kmers = 0;           // distinct locus (not private) and decoy k-mers in them
    std::uint64_t relabelled = 0;              // locus k-mers it alone holds outside the locus: private when held out
    std::uint64_t unlabelled = 0;              // decoy k-mers it alone holds: unknown when held out
    std::uint64_t copy_bases = 0;              // bases that read as a copy of the locus: not off-target
};

struct LocusIndex {
    // manifest
    std::string name, graph_path, graph_checksum, reference_path, contig;
    std::uint64_t start = 0, end = 0;
    std::string background, store_id, background_kind, coordinate_reference;
    std::uint64_t content = 0;              // the background's content_id (its genomes)
    std::uint32_t informative_max = 8;
    std::uint64_t haplotypes = 1;           // background haplotypes: S is up to informative_max times this
    std::uint64_t pad = 1500;
    double dust = 2.5;
    std::uint32_t decoy_min_hits = 20;
    std::vector<std::string> decoy_genomes;
    std::vector<std::string> calibration_genomes_given;   // --calibration-genomes: the candidates; empty: the decoy
                                                          // genomes other than the coordinate reference
    std::uint64_t calibration_max = 8;      // candidates a locus holds out, at most
    bool calibrated = false;
    std::size_t paths = 0;
    std::uint64_t decoy_windows = 0, decoy_bases = 0;
    // k-mers, sorted by code, with parallel fields
    std::vector<std::uint64_t> codes;
    std::vector<std::uint8_t> cls;
    std::vector<std::uint32_t> outside;     // saturated at 2^32 - 1
    std::vector<std::uint16_t> mult;        // largest copy number within one haplotype
    std::vector<std::uint8_t> backbone;
    std::vector<std::uint32_t> nodeset;
    std::vector<std::uint32_t> nodeset_first;   // nodesets + 1 offsets into nodeset_nodes
    std::vector<std::uint32_t> nodeset_nodes;
    std::vector<std::string> nodes;              // node names by index
    std::vector<std::uint64_t> decoy;            // sorted
    std::vector<std::uint64_t> screen;           // sorted closed syncmers among P and S
    // report
    std::vector<GenomeSpan> spans;
    std::vector<DecoyWindow> decoy_window_list;
    std::vector<DecoyWindow> calibration_window_list;   // homologous windows of calibration candidates that are not decoy
                                                        // genomes: where the anchors place the coordinate reference's
    std::vector<CalibrationGenome> calibration_genomes;  // the genomes held out in calibration, in the order chosen
    std::vector<Ceiling> ceilings;
    std::vector<Calibration> calibration;        // empty when not calibrated
    std::uint64_t dust_demoted = 0, span_demoted = 0, negative = 0;

    std::int64_t find(std::uint64_t canonical) const;
    std::string key() const;                     // what must match for a saved index to be reused
};

struct BuildOptions {
    std::vector<std::string> loci;               // locus graphs (GFA)
    std::string background;                      // store folder (its .sshash beside it is used when current)
    std::string reference_query;                 // -r: the graph path carrying the coordinates
    std::string coordinate_reference;            // the store's; default its only one
    std::uint32_t informative_max = 8;
    std::uint64_t pad = 1500;
    double dust = 2.5;
    std::uint32_t decoy_min_hits = 20;           // shared k-mer hits that make a cluster a homologous window
    std::vector<std::string> decoy_genomes;      // default: the coordinate reference
    // The candidates for the genomes held out one at a time in calibration (empty: the decoy genomes other than the
    // coordinate reference); each locus takes up to calibration_max of them, greedily by the homolog k-mers each adds.
    std::vector<std::string> calibration_genomes;
    std::size_t calibration_max = 8;
    bool calibrate = true;
    std::size_t threads = 2;
    std::function<void(const std::string&)> log;
    // Per locus, why it failed after the shared genome scans (empty: built). Given, the other loci are still
    // returned whole; without it the first failure throws.
    std::vector<std::string>* errors = nullptr;
};

// A locus's name: its file's name without .gz and .gfa, .gfa1, .fa, .fasta or .fna; its index is <name>.plx.
std::string locus_name(const std::string& path);
// The key a locus index built with these options would carry (to test a saved one before building).
std::string expected_key(const BuildOptions& o, const std::string& gfa);
std::vector<LocusIndex> build(const BuildOptions& options);

void save(const LocusIndex& index, const std::string& path);
LocusIndex load(const std::string& path);
void write_report(const LocusIndex& index, const std::string& path);

}  // namespace panvar::recruit
