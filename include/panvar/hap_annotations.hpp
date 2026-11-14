#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/gtf.hpp"

namespace panvar {

// A transcript of a haplotype's own annotation (a GFF3 or GTF in its assembly's coordinates), placed on that
// haplotype's locus path.
struct AnnotatedTranscript {
    TranscriptFeature t;          // ids and names as below; exons 0-based, half-open on the path, ascending
    std::size_t path = 0;         // index into Graph::paths
    std::string gene_record;      // the annotation's gene record: one per copy
    int copy = 0, copies = 0;     // the record's rank among the gene's records on the path, by position
    std::size_t bp_5p = 0, bp_3p = 0;   // exonic bases past the path, at each end of the transcript
    std::string detail;           // the source's own tags (CAT: source transcript and gene, class, mode)
};

// Reads `table` (haplotype<TAB>path per line, haplotype as sample#hap; a header line naming `haplotype` is
// skipped) and each annotation it names, keeping the transcripts with an exon on one of that haplotype's
// paths in `graph`. The annotation's contig names match a path's contig exactly, or once the PanSN prefix
// sample#hap# is removed, or through `alias` (a TSV whose rows list names of one contig, chromAlias-like;
// lines starting with '#' skipped). A path's interval is read as for the reference path.
//
// transcript_id is the source transcript without its version where the annotation gives one (CAT's
// source_transcript), else the transcript's transcript_id or ID; gene_name its gene's gene_name, or
// source_gene_common_name, gene or Name. In GFF3 a transcript is any feature that exons name as Parent, and
// its gene record the feature it names as Parent; in GTF they are the exons' transcript_id within their
// gene_id, and that gene_id.
std::vector<AnnotatedTranscript> load_hap_annotations(const Graph& graph, const std::string& table,
                                                      const std::string& alias, const std::string& module);

}  // namespace panvar
