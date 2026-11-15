#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/gtf.hpp"
#include "panvar/hap_annotations.hpp"

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
// to the nearest position the walk shares with the reference, and is flagged. Where a deletion leaves one copy
// of a segment the reference walks twice inside the gene, with nodes of both, the exons of the copy the pass
// does not follow are carried through a second chain of that copy's own nodes. A gene whose nodes fall short
// of the share on a walk (alleles the graph keeps on nodes of their own) is placed by sequence instead: its
// reference stretch aligned inside the window its best chain points to, accepted on the marker's coverage
// and identity.

struct TranscriptsOptions {
    std::string gfa_path;
    std::string reference_path;
    std::string gtf_path;           // reference annotation to project; empty: annotated transcripts only
    std::string hap_annotations;    // table of haplotype-specific annotations (hap_annotations.hpp)
    std::string contig_alias;
    std::string out_prefix;
    double min_exon_share = 0.8;
    double min_identity = 0.93;     // of the marker's exonic bases, for a copy placed by sequence
    std::size_t threads = 0;
    bool quiet = false;
};

// Least to most severe; a record takes the most severe that applies.
enum class HstStatus { Ok, SpliceSiteChanged, BoundaryMoved, ClippedAtEdge, IntronLost, NotAsAnnotated, MissingExon, PartialCopy, NoCopy };
constexpr int kHstStatuses = 9;
const char* hst_status_name(HstStatus s);

// Where an HST comes from: projected from the reference annotation, cut from the haplotype's own annotation,
// or both with the same sequence.
enum HstSource { kProjected = 1, kAnnotated = 2 };

// One transcript on one haplotype copy.
struct HstRecord {
    int source = kProjected;        // HstSource bits
    std::size_t transcript = 0;     // index into the projected transcripts, or (annotated only) the annotated ones
    int annotated_copies = -1;      // the gene's copies in the haplotype's annotation, where it has one
    std::string record;             // the annotation's own transcript id(s), ','-joined
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
    // Per transcript, what lies past the reference path's ends: a transcript reaching past them is projected
    // from its exons on the path, cut at its ends.
    struct EdgeClip {
        std::size_t exons = 0;        // exons in the annotation
        std::size_t first_exon = 0;   // the first one on the path, in genomic order (0-based)
        std::size_t bp_5p = 0, bp_3p = 0;   // annotated exonic bases past the path, at each end of the transcript
    };
    std::vector<EdgeClip> clips;
    std::vector<TranscriptFeature> outside;        // transcripts with no exon on the reference path
    std::vector<std::string> reference_sequences;  // per transcript, spliced from the reference
    std::vector<std::string> haplotypes;           // path names, Graph::paths order
    std::vector<AnnotatedTranscript> annotated;    // the haplotypes' own transcripts
    std::vector<char> has_annotation;              // per path: its haplotype has an annotation
    std::vector<std::string> not_given_back;       // transcripts the reference walk does not give back as annotated:
                                                   // their records on it have status not_as_annotated
};

// Projects every transcript of the GTF with an exon on the reference path onto every path of `graph`, and
// adds the transcripts of each haplotype's own annotation, cut from its walk; an annotated transcript with
// the projected one's id, copy and sequence is one record of both sources. Records reach `sink` haplotype
// by haplotype, in path order: projected ones in gene, transcript, copy order, then the annotated ones left.
class HstSink {
public:
    virtual ~HstSink() = default;
    virtual void add(const TranscriptsResult& result, const std::vector<HstRecord>& records) = 0;
};

TranscriptsResult project_transcripts(const Graph& graph, const TranscriptsOptions& options, HstSink& sink);

}  // namespace panvar
