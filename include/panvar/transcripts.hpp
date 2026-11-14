#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/gtf.hpp"

namespace panvar {

// A reference annotation projected onto every haplotype walk of a locus graph: one haplotype-specific
// transcript (HST) per transcript, haplotype and copy of its gene.
//
// A copy is one pass of a walk through the gene: a colinear chain of the nodes the walk shares with the
// gene's stretch of the reference walk, in either orientation, covering at least `min_exon_share` of the
// gene's marker exons (the coding exons of its canonical transcript, else of its longest). Genes apart on
// the reference whose passes take the same stretch of a walk are paralogs folded onto the same nodes; that
// stretch is the copy of the gene sharing more of its marker exons, and a shorter pass lying inside
// another gene's copy that matches it no better is part of that copy. Each exon boundary is carried to the
// walk through the node holding it; a boundary inside a variant (the walk misses that reference node) moves
// to the nearest position the walk shares with the reference, and is flagged.

struct TranscriptsOptions {
    std::string gfa_path;
    std::string reference_path;
    std::string gtf_path;
    std::string out_prefix;
    double min_exon_share = 0.8;
    std::size_t threads = 0;
    bool quiet = false;
};

enum class HstStatus { Ok, SpliceSiteChanged, BoundaryMoved, MissingExon, NoCopy };
const char* hst_status_name(HstStatus s);

// One transcript on one haplotype copy.
struct HstRecord {
    std::size_t transcript = 0;     // index into the projected transcripts
    std::size_t haplotype = 0;      // index into Graph::paths
    int copy = 0;                   // 1-based along the walk; 0 for no copy
    int copies = 0;                 // passes of the gene on this walk
    HstStatus status = HstStatus::NoCopy;
    bool reverse = false;           // the walk passes the gene against the reference's orientation
    std::size_t path_start = 0;     // 0-based, half-open, along the walk as written
    std::size_t path_end = 0;
    std::size_t exons = 0;          // exons projected
    int moved = 0;                  // exon boundaries moved off a variant
    int splice_changed = 0;         // splice-site dinucleotides that differ from the reference
    // Exons in transcript order: their first base along the walk as written, and their lengths.
    std::vector<std::size_t> exon_starts, exon_lengths;
    std::string detail;             // per-exon notes, ';'-joined
    std::string sequence;           // the HST, 5' to 3'; empty unless the transcript projects whole
};

struct TranscriptsResult {
    std::string reference;                         // resolved reference path name
    std::vector<TranscriptFeature> transcripts;    // projected transcripts, gene-grouped
    std::vector<TranscriptFeature> partial;        // transcripts reaching past the reference path
    std::vector<std::string> reference_sequences;  // per transcript, spliced from the reference
    std::vector<std::string> haplotypes;           // path names, Graph::paths order
};

// Projects every transcript of the GTF lying on the reference path onto every path of `graph`. Records
// reach `sink` haplotype by haplotype, in path order and gene, transcript, copy order within one.
class HstSink {
public:
    virtual ~HstSink() = default;
    virtual void add(const TranscriptsResult& result, const std::vector<HstRecord>& records) = 0;
};

TranscriptsResult project_transcripts(const Graph& graph, const TranscriptsOptions& options, HstSink& sink);

}  // namespace panvar
