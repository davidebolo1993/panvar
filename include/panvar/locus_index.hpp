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
// With the locus's haplotype-specific transcripts (panvar transcripts), their k-mers join the table: those
// the genomic walks lack (mostly across splice junctions) are flagged transcript-only. No genome holds a
// transcript, so the genomes' counts see such a k-mer only in a processed pseudogene; a paralog's transcript
// with the same junction is seen through the background's transcriptome (every annotated transcript, counted
// apart from the genomes) and through the other loci indexed in the same run. A transcript-only k-mer's
// count outside is the genomes' count (nothing subtracted), plus the transcriptome's count less its
// occurrences in the transcripts of the locus's own genes (named in the transcriptome's headers), plus one
// per other locus of the run holding it. Genomic k-mers are labelled from the genomes alone.

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
struct Calibration {
    std::string read_class;                   // short, hifi, ont
    double error_rate = 0;
    std::uint64_t read_length = 0, insert = 0;   // insert 0: single reads
    std::uint64_t held_out = 0;               // samples held out
    std::uint64_t locus = 0;                  // their fragments
    std::uint64_t flank = 0, decoy = 0;       // off-target fragments simulated
    std::uint64_t flank_screened = 0, decoy_screened = 0;
    std::array<std::uint64_t, kMaxCut + 1> locus_certain{}, flank_certain{}, decoy_certain{};
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
    std::vector<Ceiling> ceilings;
    std::vector<Calibration> calibration;        // empty when not calibrated
    std::uint64_t dust_demoted = 0, span_demoted = 0, negative = 0;
    // transcripts (empty without): their file, its checksum, whether only status=ok ones were taken or those
    // with a changed splice site too, those taken and left out (in all and by status, as text), and per k-mer
    // the transcript-only flag
    std::string transcripts, transcripts_checksum, transcripts_by_status;
    bool transcripts_ok_only = false, transcripts_splice_changed = false;
    std::uint64_t transcripts_used = 0, transcripts_skipped = 0;
    std::vector<std::uint8_t> transcript_only;
    // the background's transcriptomes when the index was built (their id, 0 without; their names), the records
    // of the locus's own genes subtracted, and the transcript-only k-mers found in other transcripts and in
    // the other loci of the run (named)
    std::uint64_t transcriptome = 0, transcriptome_records = 0, tx_in_transcriptome = 0, tx_in_other_loci = 0;
    std::string transcriptome_names, other_loci;
    // the other loci of the run it was built in (their graphs, transcripts and reference paths), hashed; 0
    // without any. Part of the key with transcripts: their k-mers decide the labels of its transcript-only
    // ones, unless they overlap it on the reference.
    std::uint64_t companions = 0;

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
    bool calibrate = true;
    // per locus (by position in `loci`): its transcripts FASTA (panvar transcripts' <prefix>.hst.fa), or empty
    std::vector<std::string> transcripts;
    bool transcripts_ok_only = false;            // take status=ok transcripts only
    bool transcripts_splice_changed = false;     // take those with a changed splice site too
    // the run's other loci that are not built now, as their saved indexes are current: what identifies each
    // (locus_token) and its index, read for the transcript-only k-mers of the loci built
    std::vector<std::string> beside_tokens, beside_indexes;
    std::size_t threads = 2;
    std::function<void(const std::string&)> log;
};

// What identifies a locus of a run to the others: its graph's checksum, its transcripts' (empty: none) and
// the reference path the query names in it, which places it on the reference.
std::string locus_token(const std::string& gfa, const std::string& transcripts, const std::string& reference_query);
// The key a locus index built with these options would carry (to test a saved one before building), for
// a locus graph with its transcripts FASTA (empty: none) and the tokens of the run's other loci.
std::string expected_key(const BuildOptions& o, const std::string& gfa, const std::string& transcripts = "",
                         const std::vector<std::string>& others = {});
std::vector<LocusIndex> build(const BuildOptions& options);

void save(const LocusIndex& index, const std::string& path);
LocusIndex load(const std::string& path);
void write_report(const LocusIndex& index, const std::string& path);

}  // namespace panvar::recruit
