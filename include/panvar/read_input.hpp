#pragma once

// A sample's reads for recruitment, as fragments (one read, or a pair of mates in sequenced orientation):
// FASTQ (plain or gzip; mates in two files or interleaved in one) and SAM, BAM or CRAM through htslib.
//
// Records: secondary, supplementary and QC-failed alignments are skipped; unmapped records and
// duplicates are kept; a reverse-strand record is turned back to the orientation it was sequenced in.
// Mates: in a file not sorted by coordinate a read waits for its mate (mates are adjacent in name-grouped
// files). A file is read as sorted by coordinate when its header says so, or says nothing and it has an
// index; positions going backwards then stop the run. In a coordinate-sorted file a read waits when its mate lies ahead within 100 kb; one whose mate
// lies further (or on another contig) is kept only when it hits a locus screen, and a read whose mate
// passed unkept is screened and, if it hits, its mate is fetched at the end, by an indexed query when the
// file has an index, else by a second pass. Either way a pair reaches the scorer when either mate hits a
// screen, as it would from FASTQ.

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace panvar::recruit {

struct Fragment {
    std::string name, seq1, seq2;
    bool paired = false;
};

struct InputOptions {
    std::string cram_reference;                               // FASTA (with .fai) for CRAM, if not resolvable otherwise
    bool interleaved = false;                                 // a single FASTQ holds mates as consecutive records
    std::size_t threads = 2;                                  // htslib decompression threads
    std::function<bool(const std::string&)> screened;         // a read hits some locus screen (coordinate-sorted input)
    std::function<void(const std::string&)> log;
};

struct InputStats {
    std::uint64_t records = 0;          // records read
    std::uint64_t secondary = 0, supplementary = 0, qc_fail = 0;   // skipped
    std::uint64_t duplicates = 0, unmapped = 0;                   // kept
    std::uint64_t fragments = 0, pairs = 0;
    std::uint64_t fetched = 0;          // mates fetched at the end (coordinate-sorted input)
    std::uint64_t unmatched = 0;        // paired reads whose mate was never found (scored alone, as single reads)
    std::uint64_t dropped = 0;          // reads set aside unscreened whose mate never came (coordinate-sorted)
    std::uint64_t unscreened = 0;       // pairs passed over because neither mate hit a screen (coordinate-sorted)
};

class FragmentSource {
public:
    virtual ~FragmentSource() = default;
    // Up to n fragments into out (resized to what was read); 0 at the end.
    virtual std::size_t next(std::vector<Fragment>& out, std::size_t n) = 0;
    virtual const InputStats& stats() const = 0;
    virtual std::string describe() const = 0;   // the file(s) and how they are read
};

// FASTQ mates in two files, read on two threads.
std::unique_ptr<FragmentSource> open_fastq_pair(const std::string& reads1, const std::string& reads2);
// One file: SAM, BAM or CRAM (by content), else FASTQ (single-end, or interleaved mates).
std::unique_ptr<FragmentSource> open_input(const std::string& path, const InputOptions& options);
// Whether htslib reads the file as SAM, BAM or CRAM.
bool is_alignment_file(const std::string& path);
// Every record of a file counted by kind (for --dry-run): `pairs` then counts paired reads and `fragments`
// single reads, as kept by the record policy; FASTQ files are single reads unless interleaved.
InputStats count_records(const std::string& path, const InputOptions& options);

}  // namespace panvar::recruit
