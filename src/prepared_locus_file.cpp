// I3-C (EXPERIMENTAL, internal): the persistent prepared-locus index v1 -- canonical codec, bounded loading with
// cumulative allocation accounting, reconstruction through the existing constructors, and atomic publication.
// Design: tests/analysis/i3c/DESIGN.md (section letters below refer to it).
#include "panvar/prepared_locus_file.hpp"

#include "panvar/genotype_repair.hpp"
#include "panvar/syncmer.hpp"
#include "prepared_locus_parts.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string_view>
#include <sys/stat.h>
#include <unistd.h>

#ifndef PANVAR_CONFIGURE_COMMIT
#define PANVAR_CONFIGURE_COMMIT "unknown"
#endif
#ifndef PANVAR_CONFIGURE_DIRTY
#define PANVAR_CONFIGURE_DIRTY "unknown"
#endif

namespace panvar {

// ---- errors, provenance -------------------------------------------------------------------------------------------------
PreparedLocusFormatError::PreparedLocusFormatError(Kind kind, const std::string& message)
    : std::runtime_error(std::string("prepared-locus index: ") + prepared_locus_error_name(kind) + ": " + message),
      kind_(kind) {}

const char* prepared_locus_error_name(PreparedLocusFormatError::Kind kind) {
    using K = PreparedLocusFormatError::Kind;
    switch (kind) {
        case K::io: return "io";
        case K::truncated: return "truncated";
        case K::bad_magic: return "bad_magic";
        case K::unsupported_schema: return "unsupported_schema";
        case K::invalid_header: return "invalid_header";
        case K::checksum_mismatch: return "checksum_mismatch";
        case K::invalid_length: return "invalid_length";
        case K::invalid_count: return "invalid_count";
        case K::overflow: return "overflow";
        case K::invalid_enum: return "invalid_enum";
        case K::invalid_section: return "invalid_section";
        case K::unknown_required_section: return "unknown_required_section";
        case K::invalid_structure: return "invalid_structure";
        case K::incompatible_preparation: return "incompatible_preparation";
        case K::limit_exceeded: return "limit_exceeded";
        case K::reconstruction_mismatch: return "reconstruction_mismatch";
        case K::exists: return "exists";
        case K::post_publication: return "post_publication";
    }
    return "unknown";
}

PreparedLocusProvenance current_prepared_locus_provenance() {
    PreparedLocusProvenance p;
    p.writer = "panvar";
    p.configure_commit = PANVAR_CONFIGURE_COMMIT;
    p.configure_dirty = PANVAR_CONFIGURE_DIRTY;
#ifdef __VERSION__
    p.compiler = __VERSION__;
#else
    p.compiler = "unknown";
#endif
#ifdef NDEBUG
    p.build_type = "ndebug";
#else
    p.build_type = "assertions";
#endif
    return p;
}

namespace {

using Kind = PreparedLocusFormatError::Kind;
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

[[noreturn]] void fail(Kind kind, const std::string& message) { throw PreparedLocusFormatError(kind, message); }

std::atomic<std::uint64_t>& load_counter() {
    static std::atomic<std::uint64_t> count{0};
    return count;
}

constexpr char kMagic[8] = {'P', 'A', 'N', 'V', 'A', 'R', 'P', 'L'};
constexpr std::size_t kHeaderBytes = 64;
constexpr std::size_t kChecksumAt = 24, kChecksumBytes = 32;
enum SectionId : std::uint32_t {
    kIdentity = 1, kCompatibility = 2, kProvenance = 3, kStructuralInput = 4, kLocusReport = 5, kNodeSequences = 6,
    kDerivedDigests = 7
};
constexpr std::uint32_t kRequiredFlag = 1;
constexpr std::uint32_t kFirstOptionalId = 16;
constexpr std::uint32_t kPanelMarkerTable = 16;   // I3-D: optional; tests/analysis/i3d/DESIGN.md

std::atomic<std::uint64_t>& precount_counter() {
    static std::atomic<std::uint64_t> count{0};
    return count;
}

std::string describe_provenance(const PreparedLocusProvenance& p) {
    return "writer=" + p.writer + ";configure_commit=" + p.configure_commit + ";configure_dirty=" + p.configure_dirty +
           ";compiler=" + p.compiler + ";build_type=" + p.build_type;
}

// ---- canonical primitives: one writer template for the file and for the content digests --------------------------------
struct StringSink {
    std::string out;
    void bytes(const void* p, std::size_t n) { out.append(static_cast<const char*>(p), n); }
};
struct HashSink {
    Sha256Stream hash;
    void bytes(const void* p, std::size_t n) { hash.update(p, n); }
};

template <typename Sink>
struct Writer {
    Sink& sink;
    void u8(std::uint8_t v) { sink.bytes(&v, 1); }
    void u32(std::uint32_t v) {
        unsigned char b[4];
        for (int i = 0; i < 4; ++i) b[i] = static_cast<unsigned char>(v >> (8 * i));
        sink.bytes(b, 4);
    }
    void u64(std::uint64_t v) {
        unsigned char b[8];
        for (int i = 0; i < 8; ++i) b[i] = static_cast<unsigned char>(v >> (8 * i));
        sink.bytes(b, 8);
    }
    void boolean(bool v) { u8(v ? 1 : 0); }
    void str(const std::string& s) {
        u64(s.size());
        sink.bytes(s.data(), s.size());
    }
    void strs(const std::vector<std::string>& v) {
        u64(v.size());
        for (const std::string& s : v) str(s);
    }
};

std::uint8_t enum_code(GenotypeBlockKind v) { return v == GenotypeBlockKind::Variant ? 0 : 1; }
std::uint8_t enum_code(ChainOrientation v) { return v == ChainOrientation::Forward ? 0 : 1; }
std::uint8_t enum_code(TemplateFrameStatus v) { return v == TemplateFrameStatus::Complete ? 0 : 1; }
std::uint8_t enum_code(BlockCoverage v) {
    switch (v) {
        case BlockCoverage::Traversing: return 0;
        case BlockCoverage::EmptyInterior: return 1;
        case BlockCoverage::BridgedBypass: return 2;
        case BlockCoverage::Uncovered: return 3;
    }
    return 255;
}

// ---- sections -----------------------------------------------------------------------------------------------------------
template <typename S>
void write_identity(Writer<S>& w, const PreparedLocusRecord& r) {
    w.str(r.key.gfa_sha256);
    w.str(r.key.reference_path);
    w.str(r.key.bubbles);
    w.u64(r.key.min_variant_bp);
    w.str(r.key.curation);
    w.u64(r.key.kmer_size);
    w.u64(r.key.syncmer_s);
    w.str(r.resolved_reference);
}

template <typename S>
void write_input(Writer<S>& w, const StructuralLocusInput& in) {
    w.str(in.locus_id);
    w.u8(enum_code(in.chain_orientation));
    w.u64(in.blocks.size());
    for (const GenotypeBlockInput& b : in.blocks) {
        w.str(b.block_id);
        w.u8(enum_code(b.kind));
        w.u64(b.alleles.size());
        for (const BlockAlleleInput& a : b.alleles) {
            w.str(a.source_id);
            w.str(a.sequence);
        }
    }
    w.strs(in.invariant_segments);
    w.u64(in.templates.size());
    for (const PanelTemplateInput& t : in.templates) {
        w.str(t.name);
        w.u8(enum_code(t.frame_status));
        w.strs(t.allele_source_ids);
        w.u8(enum_code(t.source_orientation));
        w.str(t.partial_reason);
    }
}

template <typename S>
void write_report(Writer<S>& w, const LocusInputReport& r) {
    w.u64(r.variant_blocks);
    w.u64(r.context_blocks);
    w.u64(r.total_source_alleles);
    w.u64(r.total_canonical_alleles);
    w.u64(r.total_curation_inactive_source_alleles);
    w.u64(r.total_curation_inactive_canonical_alleles);
    w.u64(r.forward_templates);
    w.u64(r.reverse_templates);
    w.strs(r.refused_templates);
    w.u64(r.contained_bubbles_dropped);
    w.u64(r.abutting_boundaries_trimmed);
    w.u64(r.boundary_bases_moved_to_context);
    w.strs(r.representation_dependent_blocks);
    w.u64(r.blocks.size());
    for (const LocusBlockReport& b : r.blocks) {
        w.str(b.block_id);
        w.u8(enum_code(b.kind));
        w.u64(b.bubble_id);
        w.u64(b.source_alleles);
        w.u64(b.canonical_alleles);
        w.u64(b.curation_inactive_source_alleles);
        w.u64(b.curation_inactive_canonical_alleles);
        w.u64(b.min_bp);
        w.u64(b.max_bp);
        w.u64(b.bypass_templates);
        w.str(b.source_node);
        w.str(b.sink_node);
        w.u64(b.trimmed_prefix_bp);
        w.u64(b.trimmed_suffix_bp);
        w.boolean(b.representation_dependent);
    }
    w.u64(r.non_traversing_coverage.size());
    for (const TemplateBlockCoverage& c : r.non_traversing_coverage) {
        w.str(c.template_name);
        w.str(c.block_id);
        w.u8(enum_code(c.coverage));
        w.boolean(c.uncovered_on_left);
    }
    w.u64(r.partial_templates.size());
    for (const PartialTemplateReport& p : r.partial_templates) {
        w.str(p.name);
        w.str(p.reason);
        w.strs(p.uncovered_blocks);
        w.str(p.missing_terminal_anchor);
    }
    w.u64(r.curated_exclusions.size());
    for (const auto& c : r.curated_exclusions) {
        w.str(c.template_name);
        w.str(c.reason);
    }
    w.str(r.curation_digest);
}

// Every PreparedMarkerPanel field, in order: shared by the panel-content digest and section 16 (one field writer).
template <typename S>
void write_panel_fields(Writer<S>& w, const PreparedMarkerPanel& p) {
    w.u64(p.kmer_size);
    w.u64(p.syncmer_s);
    w.u64(p.codes.size());
    for (const std::uint64_t c : p.codes) w.u64(c);
    w.u64(p.blocks.size());
    for (const BlockMarkerInventory& b : p.blocks) {
        w.str(b.block_id);
        w.u64(b.allele_count);
        w.u64(b.marker_slots.size());
        for (const std::uint32_t x : b.marker_slots) w.u32(x);
        w.u64(b.fixed_haploid_multiplicity.size());
        for (const std::uint32_t x : b.fixed_haploid_multiplicity) w.u32(x);
        w.u64(b.allele_markers.size());
        for (const auto& allele : b.allele_markers) {
            w.u64(allele.size());
            for (const auto& [slot, copies] : allele) {
                w.u32(slot);
                w.u32(copies);
            }
        }
    }
    w.u64(p.anchor_slots.size());
    for (const std::uint32_t x : p.anchor_slots) w.u32(x);
    w.u64(p.anchor_haploid_multiplicity.size());
    for (const std::uint32_t x : p.anchor_haploid_multiplicity) w.u32(x);
    w.u64(p.total_alleles);
    w.u64(p.scored_alleles);
    w.u64(p.dropped_multiblock_markers);
}

std::string section_body(std::uint32_t id, const PreparedLocusRecord& r) {
    StringSink sink;
    Writer<StringSink> w{sink};
    switch (id) {
        case kIdentity: write_identity(w, r); break;
        case kCompatibility:
            w.u32(r.preparation_algorithm_version);
            w.boolean(r.round_trip_passed);
            w.boolean(r.representation_independent);
            break;
        case kProvenance:
            w.str(r.provenance.writer);
            w.str(r.provenance.configure_commit);
            w.str(r.provenance.configure_dirty);
            w.str(r.provenance.compiler);
            w.str(r.provenance.build_type);
            break;
        case kStructuralInput: write_input(w, r.input); break;
        case kLocusReport: write_report(w, r.report); break;
        case kNodeSequences:
            w.u64(r.nodes.size());
            for (const auto& [id_text, seq] : r.nodes) {
                w.str(id_text);
                w.str(seq);
            }
            break;
        case kDerivedDigests:
            w.str(r.model_fingerprint);
            w.str(r.model_content_sha256);
            w.str(r.panel_content_sha256);
            break;
        case kPanelMarkerTable:
            w.u32(kPanelExtensionVersion);
            w.u32(kMarkerPanelAlgorithmVersion);
            w.str(r.model_content_sha256);   // dependency: the model the panel was built from
            write_panel_fields(w, *r.panel);  // the panel's own kmer_size and syncmer_s lead (dependencies too)
            break;
        default: fail(Kind::invalid_section, "no encoder for section " + std::to_string(id));
    }
    return std::move(sink.out);
}

std::array<unsigned char, 32> hex_to_bytes(const std::string& hex) {
    std::array<unsigned char, 32> out{};
    auto nibble = [](char c) -> int { return c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : -1; };
    if (hex.size() != 64) fail(Kind::io, "SHA-256 hex has the wrong length");
    for (std::size_t i = 0; i < 32; ++i) {
        const int hi = nibble(hex[2 * i]), lo = nibble(hex[2 * i + 1]);
        if (hi < 0 || lo < 0) fail(Kind::io, "SHA-256 hex is malformed");
        out[i] = static_cast<unsigned char>(hi * 16 + lo);
    }
    return out;
}

std::string bytes_to_hex(const unsigned char* p, std::size_t n) {
    static const char* digits = "0123456789abcdef";
    std::string out(2 * n, '0');
    for (std::size_t i = 0; i < n; ++i) {
        out[2 * i] = digits[p[i] >> 4];
        out[2 * i + 1] = digits[p[i] & 15];
    }
    return out;
}

void put_u32(std::string& s, std::size_t at, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) s[at + i] = static_cast<char>((v >> (8 * i)) & 0xff);
}
void put_u64(std::string& s, std::size_t at, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) s[at + i] = static_cast<char>((v >> (8 * i)) & 0xff);
}
std::uint32_t get_u32(const std::string& s, std::size_t at) {
    std::uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(static_cast<unsigned char>(s[at + i])) << (8 * i);
    return v;
}
std::uint64_t get_u64(const std::string& s, std::size_t at) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(static_cast<unsigned char>(s[at + i])) << (8 * i);
    return v;
}

// SHA-256 of the complete header with its checksum field zeroed, then the payload (DESIGN A.1).
std::string file_checksum_hex(const std::string& file) {
    std::string header = file.substr(0, kHeaderBytes);
    std::fill(header.begin() + kChecksumAt, header.begin() + kChecksumAt + kChecksumBytes, '\0');
    Sha256Stream hash;
    hash.update(header.data(), header.size());
    hash.update(file.data() + kHeaderBytes, file.size() - kHeaderBytes);
    return hash.hex();
}

// ---- cumulative allocation accounting (DESIGN D) ------------------------------------------------------------------------
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (a > std::numeric_limits<std::uint64_t>::max() - b) fail(Kind::overflow, "accounting addition overflows");
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a) fail(Kind::overflow, "accounting product overflows");
    return a * b;
}
std::uint64_t round_up(std::uint64_t n, std::uint64_t q) { return mul((add(n, q - 1)) / q, q); }
// The modelled cost of one heap block: macOS libmalloc quanta plus 16 bytes of header slack (also covers glibc's).
std::uint64_t alloc(std::uint64_t n) {
    if (n == 0) return 0;
    const std::uint64_t q = n <= 1008 ? round_up(n, 16) : n <= (128u << 10) ? round_up(n, 512) : round_up(n, 16u << 10);
    return add(q, 16);
}
constexpr std::uint64_t kAllocSlack = (16u << 10) + 16;   // alloc(n) <= n + kAllocSlack for every n
std::uint64_t str_heap(std::uint64_t length) { return length > 15 ? alloc(add(length, 1)) : 0; }

class Accounting {
public:
    explicit Accounting(std::uint64_t limit) : limit_(limit) {}
    void charge(std::uint64_t bytes, const char* what) {
        const std::uint64_t next = add(used_, bytes);
        if (next > limit_) {
            fail(Kind::limit_exceeded, std::string("accounted storage would reach ") + std::to_string(next) +
                                           " bytes (limit " + std::to_string(limit_) + ") at " + what);
        }
        used_ = next;
        if (tally_ != nullptr) *tally_ = add(*tally_, bytes);
        if (phase_peaks_.empty() || phase_peaks_.back().first != what) phase_peaks_.emplace_back(what, used_);
        phase_peaks_.back().second = std::max(phase_peaks_.back().second, used_);
        if (used_ > peak_) {
            peak_ = used_;
            peak_phase_ = what;
        }
    }
    void release(std::uint64_t bytes) { used_ = bytes > used_ ? 0 : used_ - bytes; }
    std::uint64_t used() const { return used_; }
    std::uint64_t peak() const { return peak_; }
    std::uint64_t limit() const { return limit_; }
    const std::string& peak_phase() const { return peak_phase_; }
    void set_tally(std::uint64_t* tally) { tally_ = tally; }
    // consecutive charges under one label form one entry: (label, highest accounted total while charging under it)
    const std::vector<std::pair<std::string, std::uint64_t>>& phase_peaks() const { return phase_peaks_; }
private:
    std::uint64_t limit_, used_ = 0, peak_ = 0;
    std::vector<std::pair<std::string, std::uint64_t>> phase_peaks_;
    std::uint64_t* tally_ = nullptr;
    std::string peak_phase_;
};

// A bounded cursor over one section body. Every count and length is validated against the remaining bytes and the
// element limit, and charged, before anything is allocated from it.
class Cursor {
public:
    Cursor(const char* begin, std::size_t size, const PreparedLocusLoadLimits& limits, Accounting& acc, const char* section)
        : p_(begin), end_(begin + size), limits_(limits), acc_(acc), section_(section) {}
    std::size_t remaining() const { return static_cast<std::size_t>(end_ - p_); }
    void need(std::size_t n) {
        if (remaining() < n) fail(Kind::invalid_length, std::string(section_) + ": field runs past the section end");
    }
    std::uint8_t u8() {
        need(1);
        return static_cast<std::uint8_t>(*p_++);
    }
    std::uint32_t u32() {
        need(4);
        std::uint32_t v = 0;
        for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(static_cast<unsigned char>(p_[i])) << (8 * i);
        p_ += 4;
        return v;
    }
    std::uint64_t u64() {
        need(8);
        std::uint64_t v = 0;
        for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(static_cast<unsigned char>(p_[i])) << (8 * i);
        p_ += 8;
        return v;
    }
    bool boolean() {
        const std::uint8_t v = u8();
        if (v > 1) fail(Kind::invalid_enum, std::string(section_) + ": a bool is not 0 or 1");
        return v == 1;
    }
    std::uint8_t enumeration(std::uint8_t max_value, const char* what) {
        const std::uint8_t v = u8();
        if (v > max_value) fail(Kind::invalid_enum, std::string(section_) + ": " + what + " code " + std::to_string(v));
        return v;
    }
    std::size_t size_value() {   // a size_t report field
        const std::uint64_t v = u64();
        if (v > std::numeric_limits<std::size_t>::max()) fail(Kind::overflow, std::string(section_) + ": value exceeds size_t");
        return static_cast<std::size_t>(v);
    }
    // A count of elements each at least `min_bytes` encoded and `object_bytes` in memory; reserves are charged here.
    std::size_t count(std::uint64_t min_bytes, std::uint64_t object_bytes, const char* what) {
        const std::uint64_t c = u64();
        if (c > limits_.max_elements) {
            fail(Kind::limit_exceeded, std::string(section_) + ": " + what + " count " + std::to_string(c) +
                                           " exceeds max_elements");
        }
        if (c != 0 && min_bytes > std::numeric_limits<std::uint64_t>::max() / c) {
            fail(Kind::overflow, std::string(section_) + ": " + what + " count times element size overflows");
        }
        if (c * min_bytes > remaining()) {
            fail(Kind::invalid_count, std::string(section_) + ": " + what + " count " + std::to_string(c) +
                                          " cannot fit in the remaining bytes");
        }
        acc_.charge(alloc(mul(c, object_bytes)), what);
        return static_cast<std::size_t>(c);
    }
    void str(std::string& out, const char* what) {
        const std::uint64_t n = u64();
        if (n > remaining()) {
            fail(Kind::invalid_length, std::string(section_) + ": " + what + " length " + std::to_string(n) +
                                           " runs past the section end");
        }
        acc_.charge(str_heap(n), what);
        out.assign(p_, static_cast<std::size_t>(n));
        p_ += n;
    }
    void strs(std::vector<std::string>& out, const char* what) {
        const std::size_t c = count(8, sizeof(std::string), what);
        out.resize(c);
        for (std::string& s : out) str(s, what);
    }
    void finish() {
        if (p_ != end_) fail(Kind::invalid_length, std::string(section_) + ": trailing bytes inside the section");
    }
private:
    const char* p_;
    const char* end_;
    const PreparedLocusLoadLimits& limits_;
    Accounting& acc_;
    const char* section_;
};

GenotypeBlockKind kind_from(std::uint8_t v) { return v == 0 ? GenotypeBlockKind::Variant : GenotypeBlockKind::Context; }
ChainOrientation orientation_from(std::uint8_t v) { return v == 0 ? ChainOrientation::Forward : ChainOrientation::Reverse; }
TemplateFrameStatus frame_from(std::uint8_t v) { return v == 0 ? TemplateFrameStatus::Complete : TemplateFrameStatus::Partial; }
BlockCoverage coverage_from(std::uint8_t v) {
    return v == 0 ? BlockCoverage::Traversing : v == 1 ? BlockCoverage::EmptyInterior
         : v == 2 ? BlockCoverage::BridgedBypass : BlockCoverage::Uncovered;
}

void read_identity(Cursor& c, PreparedLocusRecord& r) {
    c.str(r.key.gfa_sha256, "gfa_sha256");
    c.str(r.key.reference_path, "reference_path");
    c.str(r.key.bubbles, "bubbles");
    r.key.min_variant_bp = c.size_value();
    c.str(r.key.curation, "curation");
    r.key.kmer_size = c.size_value();
    r.key.syncmer_s = c.size_value();
    c.str(r.resolved_reference, "resolved_reference");
}

void read_input(Cursor& c, StructuralLocusInput& in) {
    c.str(in.locus_id, "locus_id");
    in.chain_orientation = orientation_from(c.enumeration(1, "chain orientation"));
    in.blocks.resize(c.count(8 + 1 + 8, sizeof(GenotypeBlockInput), "blocks"));
    for (GenotypeBlockInput& b : in.blocks) {
        c.str(b.block_id, "block_id");
        b.kind = kind_from(c.enumeration(1, "block kind"));
        b.alleles.resize(c.count(16, sizeof(BlockAlleleInput), "alleles"));
        for (BlockAlleleInput& a : b.alleles) {
            c.str(a.source_id, "source_id");
            c.str(a.sequence, "allele sequence");
        }
    }
    c.strs(in.invariant_segments, "invariant segments");
    in.templates.resize(c.count(8 + 1 + 8 + 1 + 8, sizeof(PanelTemplateInput), "templates"));
    for (PanelTemplateInput& t : in.templates) {
        c.str(t.name, "template name");
        t.frame_status = frame_from(c.enumeration(1, "frame status"));
        c.strs(t.allele_source_ids, "template source ids");
        t.source_orientation = orientation_from(c.enumeration(1, "template orientation"));
        c.str(t.partial_reason, "partial reason");
    }
}

void read_report(Cursor& c, LocusInputReport& r) {
    r.variant_blocks = c.size_value();
    r.context_blocks = c.size_value();
    r.total_source_alleles = c.size_value();
    r.total_canonical_alleles = c.size_value();
    r.total_curation_inactive_source_alleles = c.size_value();
    r.total_curation_inactive_canonical_alleles = c.size_value();
    r.forward_templates = c.size_value();
    r.reverse_templates = c.size_value();
    c.strs(r.refused_templates, "refused templates");
    r.contained_bubbles_dropped = c.size_value();
    r.abutting_boundaries_trimmed = c.size_value();
    r.boundary_bases_moved_to_context = c.size_value();
    c.strs(r.representation_dependent_blocks, "representation-dependent blocks");
    r.blocks.resize(c.count(8 + 1 + 8 * 8 + 8 + 8 + 8 * 2 + 1, sizeof(LocusBlockReport), "block reports"));
    for (LocusBlockReport& b : r.blocks) {
        c.str(b.block_id, "block report id");
        b.kind = kind_from(c.enumeration(1, "block report kind"));
        b.bubble_id = c.size_value();
        b.source_alleles = c.size_value();
        b.canonical_alleles = c.size_value();
        b.curation_inactive_source_alleles = c.size_value();
        b.curation_inactive_canonical_alleles = c.size_value();
        b.min_bp = c.size_value();
        b.max_bp = c.size_value();
        b.bypass_templates = c.size_value();
        c.str(b.source_node, "source node");
        c.str(b.sink_node, "sink node");
        b.trimmed_prefix_bp = c.size_value();
        b.trimmed_suffix_bp = c.size_value();
        b.representation_dependent = c.boolean();
    }
    r.non_traversing_coverage.resize(c.count(8 + 8 + 1 + 1, sizeof(TemplateBlockCoverage), "coverage rows"));
    for (TemplateBlockCoverage& row : r.non_traversing_coverage) {
        c.str(row.template_name, "coverage template");
        c.str(row.block_id, "coverage block");
        row.coverage = coverage_from(c.enumeration(3, "block coverage"));
        row.uncovered_on_left = c.boolean();
    }
    r.partial_templates.resize(c.count(8 * 4, sizeof(PartialTemplateReport), "partial templates"));
    for (PartialTemplateReport& p : r.partial_templates) {
        c.str(p.name, "partial name");
        c.str(p.reason, "partial reason");
        c.strs(p.uncovered_blocks, "uncovered blocks");
        c.str(p.missing_terminal_anchor, "missing terminal anchor");
    }
    r.curated_exclusions.resize(c.count(16, sizeof(LocusInputReport::CuratedExclusion), "curated exclusions"));
    for (auto& e : r.curated_exclusions) {
        c.str(e.template_name, "curated name");
        c.str(e.reason, "curated reason");
    }
    c.str(r.curation_digest, "curation digest");
}

// I3-D: section 16's panel fields (the order of write_panel_fields), decoded into the value the locus will own. Every
// vector is reserved exactly from its validated count, which Cursor::count charged before the allocation.
void read_panel(Cursor& c, PreparedMarkerPanel& p) {
    p.kmer_size = c.size_value();
    p.syncmer_s = c.size_value();
    p.codes.resize(c.count(8, 8, "panel codes"));
    for (std::uint64_t& code : p.codes) code = c.u64();
    p.blocks.resize(c.count(8 + 8 + 8 + 8 + 8, sizeof(BlockMarkerInventory), "panel blocks"));
    for (BlockMarkerInventory& b : p.blocks) {
        c.str(b.block_id, "panel block id");
        b.allele_count = c.size_value();
        b.marker_slots.resize(c.count(4, 4, "panel marker slots"));
        for (std::uint32_t& x : b.marker_slots) x = c.u32();
        b.fixed_haploid_multiplicity.resize(c.count(4, 4, "panel fixed multiplicities"));
        for (std::uint32_t& x : b.fixed_haploid_multiplicity) x = c.u32();
        b.allele_markers.resize(c.count(8, sizeof(std::vector<std::pair<std::uint32_t, std::uint32_t>>), "panel alleles"));
        for (auto& allele : b.allele_markers) {
            allele.resize(c.count(8, 8, "panel allele markers"));
            for (auto& [slot, copies] : allele) {
                slot = c.u32();
                copies = c.u32();
            }
        }
    }
    p.anchor_slots.resize(c.count(4, 4, "panel anchor slots"));
    for (std::uint32_t& x : p.anchor_slots) x = c.u32();
    p.anchor_haploid_multiplicity.resize(c.count(4, 4, "panel anchor multiplicities"));
    for (std::uint32_t& x : p.anchor_haploid_multiplicity) x = c.u32();
    p.total_alleles = c.size_value();
    p.scored_alleles = c.size_value();
    p.dropped_multiblock_markers = c.size_value();
}

// ---- structural invariants (DESIGN C) -----------------------------------------------------------------------------------
[[noreturn]] void bad(const std::string& message) { fail(Kind::invalid_structure, message); }

// Sorted pointers to strings, charged while they exist (the checks below allocate nothing else).
class SortedNames {
public:
    SortedNames(Accounting& acc, std::size_t reserve, const char* what) : acc_(acc), charge_(alloc(mul(reserve, 8))) {
        acc_.charge(charge_, what);
        v_.reserve(reserve);
    }
    ~SortedNames() { acc_.release(charge_); }
    void add(const std::string& s) { v_.push_back(&s); }
    // Sorts; returns the first duplicate or nullptr.
    const std::string* sort_find_duplicate() {
        std::sort(v_.begin(), v_.end(), [](const std::string* a, const std::string* b) { return *a < *b; });
        for (std::size_t i = 1; i < v_.size(); ++i) {
            if (*v_[i] == *v_[i - 1]) return v_[i];
        }
        return nullptr;
    }
    bool contains(const std::string& s) const {
        return std::binary_search(v_.begin(), v_.end(), &s, [](const std::string* a, const std::string* b) { return *a < *b; });
    }
    std::size_t distinct() const {
        std::size_t n = 0;
        for (std::size_t i = 0; i < v_.size(); ++i) n += (i == 0 || *v_[i] != *v_[i - 1]) ? 1 : 0;
        return n;
    }
private:
    Accounting& acc_;
    std::uint64_t charge_;
    std::vector<const std::string*> v_;
};

bool node_exists(const std::vector<std::pair<std::string, std::string>>& nodes, std::string_view id) {
    const auto at = std::lower_bound(nodes.begin(), nodes.end(), id,
                                     [](const auto& n, std::string_view key) { return std::string_view(n.first) < key; });
    return at != nodes.end() && std::string_view(at->first) == id;
}

// The dosage census's grammar (genotype_dosage.cpp all_walks): text after "walk:" is comma-separated node+/node- steps.
void check_walk(const std::string& source, const std::vector<std::pair<std::string, std::string>>& nodes) {
    const std::size_t at = source.find("walk:");
    if (at == std::string::npos) return;
    const std::string_view steps = std::string_view(source).substr(at + 5);
    std::size_t begin = 0;
    while (begin <= steps.size()) {
        std::size_t end = steps.find(',', begin);
        if (end == std::string_view::npos) end = steps.size();
        const std::string_view step = steps.substr(begin, end - begin);
        if (step.size() < 2 || (step.back() != '+' && step.back() != '-')) bad("malformed walk step in source " + source);
        if (!node_exists(nodes, step.substr(0, step.size() - 1))) {
            bad("walk names node '" + std::string(step.substr(0, step.size() - 1)) + "' absent from NODE_SEQUENCES");
        }
        if (end == steps.size()) break;
        begin = end + 1;
    }
}

std::uint64_t reverse_complement(std::uint64_t code, std::size_t k) {
    std::uint64_t out = 0;
    for (std::size_t i = 0; i < k; ++i) {
        out = (out << 2) | (3 - (code & 3));
        code >>= 2;
    }
    return out;
}

// I3-D DESIGN 4: the persisted panel against the decoded base, before any use. Every index in range, every dimension
// equal to the model's, and the properties build_marker_panel guarantees by construction.
void check_panel(const PreparedLocusRecord& r, Accounting& acc) {
    const PreparedMarkerPanel& p = *r.panel;
    const std::size_t k = r.key.kmer_size, s = r.key.syncmer_s;
    auto panel_bad = [](const std::string& m) { fail(Kind::invalid_structure, "panel: " + m); };
    if (p.kmer_size != k || p.syncmer_s != s) panel_bad("kmer_size or syncmer_s differs from the identity");
    if (p.codes.size() > std::numeric_limits<std::uint32_t>::max()) panel_bad("more codes than uint32 slots");
    const std::uint64_t universe = k == 32 ? ~0ull : (1ull << (2 * k));
    for (std::size_t i = 0; i < p.codes.size(); ++i) {
        const std::uint64_t code = p.codes[i];
        if (i > 0 && !(p.codes[i - 1] < code)) panel_bad("codes are not strictly increasing");
        if (code >= universe) panel_bad("a code exceeds 4^k");
        if (code > reverse_complement(code, k)) panel_bad("a code is not canonical");
        if (!is_closed_syncmer(code, k, s)) panel_bad("a code is not a closed syncmer");
    }
    if (p.blocks.size() != r.input.blocks.size()) panel_bad("one inventory per model block is required");
    const std::uint64_t slot_charge = alloc(p.codes.size());
    acc.charge(slot_charge, "panel check: slot partition");
    std::vector<std::uint8_t> owner(p.codes.size(), 0);
    std::uint64_t alleles = 0;
    for (std::size_t bi = 0; bi < p.blocks.size(); ++bi) {
        const BlockMarkerInventory& b = p.blocks[bi];
        if (b.block_id != r.input.blocks[bi].block_id) panel_bad("block " + std::to_string(bi) + " id differs from the model");
        if (b.allele_count != r.report.blocks[bi].canonical_alleles) panel_bad("block " + b.block_id + ": allele count");
        if (b.allele_markers.size() != b.allele_count) panel_bad("block " + b.block_id + ": allele marker lists");
        if (b.marker_slots.size() != b.fixed_haploid_multiplicity.size()) panel_bad("block " + b.block_id + ": fixed multiplicities");
        alleles += b.allele_count;
        for (std::size_t j = 0; j < b.marker_slots.size(); ++j) {
            const std::uint32_t slot = b.marker_slots[j];
            if (slot >= p.codes.size()) panel_bad("block " + b.block_id + ": a marker slot is out of range");
            if (j > 0 && !(b.marker_slots[j - 1] < slot)) panel_bad("block " + b.block_id + ": marker slots are not increasing");
            if (owner[slot] != 0) panel_bad("a code is a marker of two blocks");
            owner[slot] = 1;
        }
        // per local marker: carriers, min and max copies (the block's own transient, charged)
        const std::size_t m = b.marker_slots.size();
        const std::uint64_t local_charge = alloc(mul(m, 12));
        acc.charge(local_charge, "panel check: block markers");
        std::vector<std::uint32_t> carriers(m, 0), lo(m, std::numeric_limits<std::uint32_t>::max()), hi(m, 0);
        for (const auto& allele : b.allele_markers) {
            for (std::size_t e = 0; e < allele.size(); ++e) {
                const auto [local, copies] = allele[e];
                if (local >= m) panel_bad("block " + b.block_id + ": an allele marker index is out of range");
                if (e > 0 && !(allele[e - 1].first < local)) panel_bad("block " + b.block_id + ": allele markers are not increasing");
                if (copies == 0) panel_bad("block " + b.block_id + ": an allele marker with zero copies");
                ++carriers[local];
                lo[local] = std::min(lo[local], copies);
                hi[local] = std::max(hi[local], copies);
            }
        }
        for (std::size_t j = 0; j < m; ++j) {
            if (carriers[j] == 0) panel_bad("block " + b.block_id + ": a variation marker no allele carries");
            if (carriers[j] == b.allele_count && lo[j] == hi[j]) {
                panel_bad("block " + b.block_id + ": a variation marker with equal copies in every allele");
            }
        }
        acc.release(local_charge);
    }
    if (p.anchor_slots.size() != p.anchor_haploid_multiplicity.size()) panel_bad("anchor slots and multiplicities differ in size");
    for (std::size_t j = 0; j < p.anchor_slots.size(); ++j) {
        const std::uint32_t slot = p.anchor_slots[j];
        if (slot >= p.codes.size()) panel_bad("an anchor slot is out of range");
        if (j > 0 && !(p.anchor_slots[j - 1] < slot)) panel_bad("anchor slots are not increasing");
        if (p.anchor_haploid_multiplicity[j] == 0) panel_bad("an anchor with zero copies");
        if (owner[slot] != 0) panel_bad("a code is both an anchor and a block marker");
        owner[slot] = 1;
    }
    for (const std::uint8_t o : owner) {
        if (o == 0) panel_bad("a code is neither an anchor nor a block marker");
    }
    acc.release(slot_charge);
    if (p.total_alleles != alleles || p.scored_alleles != alleles) panel_bad("total or scored alleles differ from the blocks");
}

void check_structure(const PreparedLocusRecord& r, Accounting& acc) {
    const StructuralLocusInput& in = r.input;
    const LocusInputReport& rep = r.report;
    // I3-D hardening: the marker settings, before any pre-count or reconstruction can see them
    if (r.key.syncmer_s < 1 || r.key.syncmer_s > r.key.kmer_size || r.key.kmer_size > 31) {
        bad("marker settings must satisfy 1 <= syncmer_s <= kmer_size <= 31 (kmer_size " + std::to_string(r.key.kmer_size) +
            ", syncmer_s " + std::to_string(r.key.syncmer_s) + ")");
    }
    if (r.preparation_algorithm_version != kPreparationAlgorithmVersion) {
        fail(Kind::incompatible_preparation, "preparation_algorithm_version " + std::to_string(r.preparation_algorithm_version) +
                                                 " (this binary: " + std::to_string(kPreparationAlgorithmVersion) + ")");
    }
    if (!r.round_trip_passed || !r.representation_independent) bad("the preparation-time checks did not both pass");
    if (r.resolved_reference != in.locus_id) bad("the resolved reference differs from the locus id");
    for (const std::string* digest : {&r.model_content_sha256, &r.panel_content_sha256}) {
        if (digest->size() != 64) bad("a content digest is not a SHA-256");
    }
    // nodes
    for (std::size_t i = 0; i < r.nodes.size(); ++i) {
        if (r.nodes[i].first.empty()) bad("an empty node id");
        if (i > 0 && !(r.nodes[i - 1].first < r.nodes[i].first)) bad("node ids are not strictly increasing");
    }
    // input
    if (in.blocks.empty()) bad("no blocks");
    if (in.invariant_segments.size() != in.blocks.size() + 1) bad("invariant segments are not blocks + 1");
    SortedNames block_ids(acc, in.blocks.size(), "invariant check: block ids");
    std::size_t largest_block = 0;
    for (const GenotypeBlockInput& b : in.blocks) {
        if (b.block_id.empty()) bad("an empty block id");
        if (b.alleles.empty()) bad("block " + b.block_id + " has no alleles");
        block_ids.add(b.block_id);
        largest_block = std::max(largest_block, b.alleles.size());
    }
    if (const std::string* d = block_ids.sort_find_duplicate()) bad("duplicate block id " + *d);
    SortedNames template_names(acc, in.templates.size(), "invariant check: template names");
    std::size_t complete_forward = 0, complete_reverse = 0;
    for (const PanelTemplateInput& t : in.templates) {
        if (t.name.empty()) bad("an empty template name");
        template_names.add(t.name);
    }
    if (const std::string* d = template_names.sort_find_duplicate()) bad("duplicate template name " + *d);
    // per block: unique sources, distinct sequences, report agreement; templates map to existing sources
    if (rep.blocks.size() != in.blocks.size()) bad("one block report per input block is required");
    std::size_t variant = 0, context = 0, total_source = 0, total_canonical = 0, inactive_source = 0, inactive_canonical = 0;
    std::vector<std::string> dependent;
    for (std::size_t bi = 0; bi < in.blocks.size(); ++bi) {
        const GenotypeBlockInput& b = in.blocks[bi];
        const LocusBlockReport& br = rep.blocks[bi];
        SortedNames sources(acc, b.alleles.size(), "invariant check: block sources");
        SortedNames sequences(acc, b.alleles.size(), "invariant check: block sequences");
        std::size_t min_bp = 0, max_bp = 0;
        for (std::size_t a = 0; a < b.alleles.size(); ++a) {
            const BlockAlleleInput& allele = b.alleles[a];
            if (allele.source_id.empty()) bad("an empty source id in block " + b.block_id);
            sources.add(allele.source_id);
            sequences.add(allele.sequence);
            check_walk(allele.source_id, r.nodes);
            min_bp = a == 0 ? allele.sequence.size() : std::min(min_bp, allele.sequence.size());
            max_bp = std::max(max_bp, allele.sequence.size());
        }
        if (const std::string* d = sources.sort_find_duplicate()) bad("duplicate source " + *d + " in block " + b.block_id);
        sequences.sort_find_duplicate();
        for (const PanelTemplateInput& t : in.templates) {
            if (t.frame_status != TemplateFrameStatus::Complete) continue;
            if (t.allele_source_ids.size() != in.blocks.size()) bad("template " + t.name + " does not map every block");
            if (!sources.contains(t.allele_source_ids[bi])) {
                bad("template " + t.name + " names a source absent from block " + b.block_id);
            }
        }
        if (br.block_id != b.block_id || br.kind != b.kind) bad("block report " + br.block_id + " differs from its block");
        if (b.kind == GenotypeBlockKind::Context && br.bubble_id != 0) bad("a context block report carries a bubble id");
        if (br.source_alleles != b.alleles.size()) bad("block " + b.block_id + ": source allele count");
        if (br.canonical_alleles != sequences.distinct()) bad("block " + b.block_id + ": canonical allele count");
        if (br.min_bp != min_bp || br.max_bp != max_bp) bad("block " + b.block_id + ": allele length range");
        for (const std::string* node : {&br.source_node, &br.sink_node}) {
            if (!node->empty() && !node_exists(r.nodes, *node)) bad("block " + b.block_id + ": boundary node absent");
        }
        (b.kind == GenotypeBlockKind::Variant ? variant : context) += 1;
        total_source += br.source_alleles;
        total_canonical += br.canonical_alleles;
        inactive_source += br.curation_inactive_source_alleles;
        inactive_canonical += br.curation_inactive_canonical_alleles;
        if (br.representation_dependent) dependent.push_back(br.block_id);
    }
    for (const PanelTemplateInput& t : in.templates) {
        if (t.frame_status == TemplateFrameStatus::Complete) {
            (t.source_orientation == ChainOrientation::Forward ? complete_forward : complete_reverse) += 1;
        }
    }
    if (rep.variant_blocks != variant || rep.context_blocks != context) bad("report block counts differ from block kinds");
    if (rep.total_source_alleles != total_source || rep.total_canonical_alleles != total_canonical ||
        rep.total_curation_inactive_source_alleles != inactive_source ||
        rep.total_curation_inactive_canonical_alleles != inactive_canonical) {
        bad("report allele totals differ from the per-block sums");
    }
    if (rep.representation_dependent_blocks != dependent) bad("representation-dependent blocks differ from the flags");
    try {
        require_representation_independent(rep);
    } catch (const std::exception& e) {
        bad(e.what());
    }
    // refused == partial (same order), and as a set the Partial input templates; curated names disjoint from templates
    if (rep.refused_templates.size() != rep.partial_templates.size()) bad("refused and partial template records differ");
    for (std::size_t i = 0; i < rep.refused_templates.size(); ++i) {
        if (rep.refused_templates[i] != rep.partial_templates[i].name) bad("refused and partial template records differ");
    }
    {
        SortedNames refused(acc, rep.refused_templates.size(), "invariant check: refused");
        for (const std::string& n : rep.refused_templates) refused.add(n);
        if (refused.sort_find_duplicate() != nullptr) bad("a refused template is recorded twice");
        std::size_t partial_inputs = 0;
        for (const PanelTemplateInput& t : in.templates) {
            if (t.frame_status != TemplateFrameStatus::Partial) continue;
            ++partial_inputs;
            if (!refused.contains(t.name)) bad("partial template " + t.name + " is not recorded as refused");
        }
        if (partial_inputs != rep.refused_templates.size()) bad("refused templates are not exactly the partial templates");
        for (const std::string& n : rep.refused_templates) {
            if (!template_names.contains(n)) bad("refused template " + n + " is not an input template");
        }
    }
    SortedNames curated(acc, rep.curated_exclusions.size(), "invariant check: curated");
    for (const auto& e : rep.curated_exclusions) {
        if (e.template_name.empty()) bad("an empty curated name");
        if (template_names.contains(e.template_name)) bad("curated template " + e.template_name + " is also an input template");
        curated.add(e.template_name);
    }
    if (const std::string* d = curated.sort_find_duplicate()) bad("curated template " + *d + " recorded twice");
    if (rep.curated_exclusions.empty() != rep.curation_digest.empty()) bad("curation digest present without exclusions, or absent with them");
    if (r.key.curation == "none" && !rep.curated_exclusions.empty()) bad("exclusions recorded without a curation manifest");
    // orientation counts are taken before curation over every oriented path (DESIGN C)
    const std::uint64_t oriented = static_cast<std::uint64_t>(rep.forward_templates) + rep.reverse_templates;
    if (oriented < complete_forward + complete_reverse || oriented > in.templates.size() + rep.curated_exclusions.size()) {
        bad("forward + reverse template counts are outside what the input and curation allow");
    }
    if (complete_forward > rep.forward_templates || complete_reverse > rep.reverse_templates) {
        bad("more complete templates of one orientation than the report counted");
    }
    for (const TemplateBlockCoverage& c : rep.non_traversing_coverage) {
        if (!template_names.contains(c.template_name) && !curated.contains(c.template_name)) {
            bad("coverage row names template " + c.template_name + ", neither an input template nor curated");
        }
        if (!block_ids.contains(c.block_id)) bad("coverage row names block " + c.block_id + ", absent");
        if (c.uncovered_on_left && c.coverage != BlockCoverage::Uncovered) bad("uncovered_on_left without Uncovered");
    }
    for (const PartialTemplateReport& p : rep.partial_templates) {
        for (const std::string& block : p.uncovered_blocks) {
            if (!block_ids.contains(block)) bad("partial template " + p.name + " names block " + block + ", absent");
        }
        const std::string& m = p.missing_terminal_anchor;
        if (!(m.empty() || m == "left" || m == "right" || m == "both")) bad("missing terminal anchor '" + m + "'");
    }
    if (r.panel) check_panel(r, acc);
}

// ---- the decode stage ----------------------------------------------------------------------------------------------------
struct DecodeTallies {
    std::uint64_t input = 0, report = 0, nodes = 0, other = 0, panel = 0;
};

PreparedLocusRecord decode_with(const std::string& file, const PreparedLocusLoadLimits& limits, Accounting& acc,
                                DecodeTallies& tallies, PreparedLocusLoadReport* report) {
    if (file.size() < kHeaderBytes) fail(Kind::truncated, "shorter than the 64-byte header");
    if (std::memcmp(file.data(), kMagic, 8) != 0) fail(Kind::bad_magic, "not a prepared-locus index");
    if (get_u32(file, 8) != kPreparedLocusSchemaVersion) {
        fail(Kind::unsupported_schema, "schema version " + std::to_string(get_u32(file, 8)));
    }
    if (get_u32(file, 12) != kHeaderBytes) fail(Kind::invalid_header, "header_bytes is not 64");
    const std::uint64_t payload = get_u64(file, 16);
    if (payload != file.size() - kHeaderBytes) {
        fail(payload > file.size() - kHeaderBytes ? Kind::truncated : Kind::invalid_header,
             "payload_bytes " + std::to_string(payload) + " but " + std::to_string(file.size() - kHeaderBytes) + " follow");
    }
    const std::uint32_t section_count = get_u32(file, 56);
    Clock::time_point mark = Clock::now();
    const std::string stored = bytes_to_hex(reinterpret_cast<const unsigned char*>(file.data()) + kChecksumAt, kChecksumBytes);
    const std::string computed = file_checksum_hex(file);
    if (stored != computed) fail(Kind::checksum_mismatch, "header and payload do not match their SHA-256");
    if (get_u32(file, 60) != 0) fail(Kind::invalid_header, "reserved header flags are set");
    if (report != nullptr) {
        report->seconds.emplace_back("load_checksum", since(mark));
        report->file_checksum = stored;
        report->file_bytes = file.size();
    }
    mark = Clock::now();
    PreparedLocusRecord r;
    std::size_t at = kHeaderBytes;
    std::uint32_t seen = 0, previous_id = 0, next_required = kIdentity;
    while (at < file.size()) {
        if (file.size() - at < 16) fail(Kind::invalid_section, "a section header runs past the payload");
        const std::uint32_t id = get_u32(file, at), flags = get_u32(file, at + 4);
        const std::uint64_t bytes = get_u64(file, at + 8);
        at += 16;
        if (bytes > file.size() - at) fail(Kind::invalid_section, "section " + std::to_string(id) + " runs past the payload");
        if ((flags & ~kRequiredFlag) != 0) fail(Kind::invalid_section, "reserved section flags set in section " + std::to_string(id));
        if (id == 0 || (seen > 0 && id <= previous_id)) fail(Kind::invalid_section, "section ids are not strictly increasing");
        previous_id = id;
        ++seen;
        const bool required = (flags & kRequiredFlag) != 0;
        if (id <= kDerivedDigests) {
            if (!required) fail(Kind::invalid_section, "v1 section " + std::to_string(id) + " is not flagged required");
            if (id != next_required) fail(Kind::invalid_section, "required section " + std::to_string(next_required) + " is missing");
            ++next_required;
            const char* names[] = {"", "IDENTITY", "COMPATIBILITY", "PROVENANCE", "STRUCTURAL_INPUT", "LOCUS_REPORT",
                                   "NODE_SEQUENCES", "DERIVED_DIGESTS"};
            Cursor c(file.data() + at, static_cast<std::size_t>(bytes), limits, acc, names[id]);
            switch (id) {
                case kIdentity: acc.set_tally(&tallies.other); read_identity(c, r); break;
                case kCompatibility:
                    acc.set_tally(&tallies.other);
                    r.preparation_algorithm_version = c.u32();
                    r.round_trip_passed = c.boolean();
                    r.representation_independent = c.boolean();
                    break;
                case kProvenance:
                    acc.set_tally(&tallies.other);
                    c.str(r.provenance.writer, "writer");
                    c.str(r.provenance.configure_commit, "configure_commit");
                    c.str(r.provenance.configure_dirty, "configure_dirty");
                    c.str(r.provenance.compiler, "compiler");
                    c.str(r.provenance.build_type, "build_type");
                    break;
                case kStructuralInput: acc.set_tally(&tallies.input); read_input(c, r.input); break;
                case kLocusReport: acc.set_tally(&tallies.report); read_report(c, r.report); break;
                case kNodeSequences: {
                    acc.set_tally(&tallies.nodes);
                    r.nodes.resize(c.count(16, sizeof(std::pair<std::string, std::string>), "nodes"));
                    for (auto& n : r.nodes) {
                        c.str(n.first, "node id");
                        c.str(n.second, "node sequence");
                    }
                    break;
                }
                case kDerivedDigests:
                    acc.set_tally(&tallies.other);
                    c.str(r.model_fingerprint, "model fingerprint");
                    c.str(r.model_content_sha256, "model content digest");
                    c.str(r.panel_content_sha256, "panel content digest");
                    break;
            }
            acc.set_tally(nullptr);
            c.finish();
        } else if (id == kPanelMarkerTable && !required && bytes >= 8 && get_u32(file, at) == kPanelExtensionVersion &&
                   get_u32(file, at + 4) == kMarkerPanelAlgorithmVersion) {
            // I3-D: understood. Malformed content refuses; it is never silently rebuilt.
            Cursor c(file.data() + at, static_cast<std::size_t>(bytes), limits, acc, "PANEL_MARKER_TABLE");
            acc.set_tally(&tallies.panel);
            (void)c.u32();
            (void)c.u32();
            std::string model_dependency;
            c.str(model_dependency, "panel model dependency");
            if (model_dependency != r.model_content_sha256) {
                fail(Kind::invalid_structure, "the panel section depends on a different model than DERIVED_DIGESTS records");
            }
            r.panel.emplace();
            read_panel(c, *r.panel);
            acc.set_tally(nullptr);
            c.finish();
        } else if (id == kPanelMarkerTable && required) {
            fail(Kind::invalid_section, "section 16 (the panel table) is optional by definition but is flagged required");
        } else if (required) {
            fail(Kind::unknown_required_section, "section " + std::to_string(id) + " is required but not understood");
        } else if (id < kFirstOptionalId) {
            fail(Kind::invalid_section, "section id " + std::to_string(id) + " is reserved");
        } else {
            // an optional section this binary does not understand (an unknown id, or section 16 at an unsupported
            // version): not used, kept verbatim so a rewrite preserves it
            acc.set_tally(&tallies.other);
            acc.charge(add(str_heap(bytes), alloc(sizeof(std::pair<std::uint32_t, std::string>) * 4)), "unknown optional section");
            acc.set_tally(nullptr);
            r.unknown_optional_sections.emplace_back(id, file.substr(at, static_cast<std::size_t>(bytes)));
        }
        at += static_cast<std::size_t>(bytes);
    }
    if (next_required != kDerivedDigests + 1) fail(Kind::invalid_section, "required section " + std::to_string(next_required) + " is missing");
    if (seen != section_count) fail(Kind::invalid_section, "section_count does not match the sections present");
    if (report != nullptr) report->seconds.emplace_back("load_decode", since(mark));
    mark = Clock::now();
    check_structure(r, acc);
    if (report != nullptr) report->seconds.emplace_back("load_invariants", since(mark));
    return r;
}

// DESIGN D phase 4: prepare_genotype_model's persistent storage and its largest per-block transient, from the input.
struct ModelProjection {
    std::uint64_t persistent = 0, transient = 0;
};
ModelProjection project_model(const StructuralLocusInput& in) {
    using MapNode = std::pair<const std::string, std::size_t>;
    using SeqNode = std::pair<const std::string, std::vector<std::string>>;
    ModelProjection p;
    const std::uint64_t B = in.blocks.size(), T = in.templates.size();
    p.persistent = alloc(mul(B, sizeof(PreparedGenotypeBlock))) + str_heap(in.locus_id.size());
    std::uint64_t worst_block = 0;
    for (const GenotypeBlockInput& b : in.blocks) {
        const std::uint64_t n = b.alleles.size();
        std::uint64_t persistent = add(str_heap(b.block_id.size()), alloc(mul(n, sizeof(CanonicalBlockAllele))));
        std::uint64_t transient = 0;
        const std::uint64_t stable_id = str_heap(b.block_id.size() + 8 + 20);
        for (const BlockAlleleInput& a : b.alleles) {
            // bounded by input alleles: stable id, sequence copy, source-id vector share (<= 96 B per source, below),
            // the source-id copy, and the allele_id_by_source node with its key copy
            persistent = add(persistent, stable_id + str_heap(a.sequence.size()) + 96 + str_heap(a.source_id.size()) +
                                             alloc(32 + sizeof(MapNode)) + str_heap(a.source_id.size()));
            // set<string> of source ids; map<string, vector<string>> keyed by sequence copies
            transient = add(transient, alloc(32 + sizeof(std::string)) + str_heap(a.source_id.size()) +
                                           alloc(32 + sizeof(SeqNode)) + str_heap(a.sequence.size()));
        }
        p.persistent = add(p.persistent, persistent);
        worst_block = std::max(worst_block, transient);
    }
    // templates: ordered pointers (transient), names, rows, orientations; validate_panel's set (transient)
    std::uint64_t names = alloc(mul(T, sizeof(std::string))) + alloc(mul(T, sizeof(std::vector<std::size_t>))) +
                          alloc(mul(T, sizeof(ChainOrientation)));
    std::uint64_t name_set = alloc(mul(T, 8));
    for (const PanelTemplateInput& t : in.templates) {
        names = add(names, str_heap(t.name.size()) + alloc(mul(B, sizeof(std::size_t))));
        name_set = add(name_set, alloc(32 + sizeof(std::string)) + str_heap(t.name.size()));
    }
    p.persistent = add(p.persistent, names);
    // invariant segments copied into the model
    p.persistent = add(p.persistent, alloc(mul(in.invariant_segments.size(), sizeof(std::string))));
    for (const std::string& seg : in.invariant_segments) p.persistent = add(p.persistent, str_heap(seg.size()));
    p.transient = add(worst_block, name_set);
    return p;
}

std::uint64_t panel_storage(const PreparedMarkerPanel& panel) {   // the built panel's actual storage (capacities)
    std::uint64_t s = alloc(mul(panel.codes.capacity(), 8)) + alloc(mul(panel.anchor_slots.capacity(), 4)) +
                      alloc(mul(panel.anchor_haploid_multiplicity.capacity(), 4)) +
                      alloc(mul(panel.blocks.capacity(), sizeof(BlockMarkerInventory)));
    for (const BlockMarkerInventory& b : panel.blocks) {
        s = add(s, str_heap(b.block_id.size()) + alloc(mul(b.marker_slots.capacity(), 4)) +
                       alloc(mul(b.fixed_haploid_multiplicity.capacity(), 4)) +
                       alloc(mul(b.allele_markers.capacity(), sizeof(std::vector<std::pair<std::uint32_t, std::uint32_t>>))));
        for (const auto& a : b.allele_markers) s = add(s, alloc(mul(a.capacity(), 8)));
    }
    return s;
}

// DESIGN D phase 7: build_marker_panel's simultaneous storage, bounded from exact syncmer occurrence counts.
struct PanelBound {
    std::uint64_t bytes = 0, occurrences = 0;
};
PanelBound bound_panel(const PreparedGenotypeModel& model, std::size_t k, std::size_t s, Accounting& acc) {
    PanelBound out;
    // the longest sequence's collect_syncmers transient, charged before any counting (static, from its length)
    std::uint64_t longest = 0;
    for (const PreparedGenotypeBlock& b : model.blocks()) {
        for (const CanonicalBlockAllele& a : b.alleles) longest = std::max<std::uint64_t>(longest, a.sequence.size());
    }
    for (const std::string& seg : model.invariant_segments()) longest = std::max<std::uint64_t>(longest, seg.size());
    const std::uint64_t positions_longest = longest >= k ? longest - k + 1 : 0;
    const std::uint64_t scan = alloc(mul(16, std::max(longest / 4 + 1, mul(2, positions_longest))));
    acc.charge(scan, "panel bound: pre-count scan");
    ++precount_counter();   // I3-D: the loaded-panel route never reaches this
    std::uint64_t O_A = 0, O_C = 0, counts = alloc(mul(model.blocks().size(), sizeof(std::vector<int>)));
    std::uint64_t markers = alloc(mul(model.blocks().size(), sizeof(BlockMarkerInventory)));
    std::uint64_t seq_transient = 0, max_block_occ = 0, max_block_alleles = 0;
    auto per_sequence = [&](std::uint64_t length, std::uint64_t o) {
        seq_transient = std::max(seq_transient, alloc(mul(16, std::max(length / 4 + 1, mul(2, o)))) + alloc(mul(8, o)) +
                                                    alloc(mul(32, o)));
    };
    for (const PreparedGenotypeBlock& b : model.blocks()) {
        const std::uint64_t n = b.alleles.size();
        counts = add(counts, alloc(mul(n, sizeof(std::vector<std::pair<std::uint64_t, std::uint32_t>>))));
        markers = add(markers, str_heap(b.block_id.size()) +
                                   alloc(mul(n, sizeof(std::vector<std::pair<std::uint32_t, std::uint32_t>>))));
        std::uint64_t block_occ = 0;
        for (const CanonicalBlockAllele& a : b.alleles) {
            const std::uint64_t o = collect_syncmers(a.sequence, k, s).size();
            O_A = add(O_A, o);
            block_occ = add(block_occ, o);
            counts = add(counts, alloc(mul(32, o)));     // CountedMarkers: <= 2 * distinct <= 2 * o entries of 16 B
            markers = add(markers, alloc(mul(16, o)));   // allele marker entries: <= 2 * o entries of 8 B
            per_sequence(a.sequence.size(), o);
        }
        max_block_occ = std::max(max_block_occ, block_occ);
        max_block_alleles = std::max(max_block_alleles, n);
    }
    for (const std::string& seg : model.invariant_segments()) {
        const std::uint64_t o = collect_syncmers(seg, k, s).size();
        O_C = add(O_C, o);
        per_sequence(seg.size(), o);
    }
    acc.release(scan);
    // Joint containers, each with capacity < 2x its entries: with V variation events and Fv fixed block entries,
    // V + Fv <= O_A; retained R <= V; aggregated fixed F' <= Fv + O_C. The bound takes the worse extreme split plus the
    // per-allocation rounding slack (alloc(n) <= n + slack), which covers any interior split.
    auto end_state = [&](std::uint64_t V, std::uint64_t Fv) {
        const std::uint64_t R = V, F = add(Fv, O_C);
        return alloc(mul(32, V)) + alloc(mul(16, V)) + alloc(mul(32, R)) + alloc(mul(32, F)) + alloc(mul(8, F)) +
               alloc(mul(8, F)) + alloc(mul(8, add(F, R))) + alloc(mul(8, R)) + alloc(mul(8, R)) + alloc(mul(8, F)) +
               alloc(mul(8, F));
    };
    auto aggregate_state = [&](std::uint64_t V, std::uint64_t Fv) {
        const std::uint64_t F = add(Fv, O_C);
        return alloc(mul(32, V)) + alloc(mul(32, add(Fv, O_C))) + alloc(mul(32, F));
    };
    const std::uint64_t blocks = model.blocks().size();
    const std::uint64_t J4 = add(std::max(end_state(O_A, 0), end_state(0, O_A)), mul(add(11, mul(2, blocks)), kAllocSlack));
    const std::uint64_t J3 = add(std::max(aggregate_state(O_A, 0), aggregate_state(0, O_A)), mul(3, kAllocSlack));
    const std::uint64_t transient = add(add(seq_transient, alloc(mul(16, max_block_occ))), alloc(mul(48, max_block_alleles)));
    out.bytes = add(std::max(add(counts, J3), add(add(counts, markers), J4)), transient);
    out.occurrences = add(O_A, O_C);
    return out;
}

void model_digest_fields(Writer<HashSink>& w, const PreparedGenotypeModel& m) {
    w.str("panvar.prepared_locus.model_content.v1");
    w.str(m.locus_id());
    w.u8(enum_code(m.chain_orientation()));
    w.u8(enum_code(m.frame_status()));
    w.u64(m.blocks().size());
    for (const PreparedGenotypeBlock& b : m.blocks()) {
        w.str(b.block_id);
        w.u8(enum_code(b.kind));
        w.u64(b.alleles.size());
        for (const CanonicalBlockAllele& a : b.alleles) {
            w.u64(a.id);
            w.str(a.stable_id);
            w.str(a.sequence);
            w.strs(a.source_ids);
            w.u64(a.panel_carrier_count);
        }
        w.u64(b.allele_id_by_source.size());
        for (const auto& [source, id] : b.allele_id_by_source) {
            w.str(source);
            w.u64(id);
        }
    }
    w.strs(m.invariant_segments());
    w.strs(m.panel().template_names);
    w.u64(m.panel().template_alleles.size());
    for (const auto& row : m.panel().template_alleles) {
        w.u64(row.size());
        for (const std::size_t id : row) w.u64(id);
    }
    w.u64(m.template_source_orientations().size());
    for (const ChainOrientation o : m.template_source_orientations()) w.u8(enum_code(o));
}

}  // namespace

std::uint64_t prepared_locus_loads() { return load_counter().load(); }

std::string prepared_locus_model_content_sha256(const PreparedGenotypeModel& model) {
    HashSink sink;
    Writer<HashSink> w{sink};
    model_digest_fields(w, model);
    return sink.hash.hex();
}

std::string prepared_locus_panel_content_sha256(const PreparedMarkerPanel& p) {
    HashSink sink;
    Writer<HashSink> w{sink};
    w.str("panvar.prepared_locus.panel_content.v1");
    write_panel_fields(w, p);
    return sink.hash.hex();
}

std::uint64_t prepared_locus_panel_precounts() { return precount_counter().load(); }

std::string prepared_marker_panel_difference(const PreparedMarkerPanel& a, const PreparedMarkerPanel& b) {
    if (a.kmer_size != b.kmer_size) return "kmer_size";
    if (a.syncmer_s != b.syncmer_s) return "syncmer_s";
    if (a.codes != b.codes) return "codes";
    if (a.blocks.size() != b.blocks.size()) return "blocks (count)";
    for (std::size_t i = 0; i < a.blocks.size(); ++i) {
        const BlockMarkerInventory &x = a.blocks[i], &y = b.blocks[i];
        const std::string at = "blocks[" + std::to_string(i) + "].";
        if (x.block_id != y.block_id) return at + "block_id";
        if (x.allele_count != y.allele_count) return at + "allele_count";
        if (x.marker_slots != y.marker_slots) return at + "marker_slots";
        if (x.fixed_haploid_multiplicity != y.fixed_haploid_multiplicity) return at + "fixed_haploid_multiplicity";
        if (x.allele_markers != y.allele_markers) return at + "allele_markers";
    }
    if (a.anchor_slots != b.anchor_slots) return "anchor_slots";
    if (a.anchor_haploid_multiplicity != b.anchor_haploid_multiplicity) return "anchor_haploid_multiplicity";
    if (a.total_alleles != b.total_alleles) return "total_alleles";
    if (a.scored_alleles != b.scored_alleles) return "scored_alleles";
    if (a.dropped_multiblock_markers != b.dropped_multiblock_markers) return "dropped_multiblock_markers";
    return "";
}

std::string encode_prepared_locus(const PreparedLocusRecord& record) {
    const std::uint32_t ids[] = {kIdentity, kCompatibility, kProvenance, kStructuralInput, kLocusReport, kNodeSequences,
                                 kDerivedDigests};
    std::string file(kHeaderBytes, '\0');
    std::uint32_t sections = 0;
    const auto append = [&](std::uint32_t id, std::uint32_t flags, const std::string& body) {
        std::string head(16, '\0');
        put_u32(head, 0, id);
        put_u32(head, 4, flags);
        put_u64(head, 8, body.size());
        file += head;
        file += body;
        ++sections;
    };
    for (const std::uint32_t id : ids) append(id, kRequiredFlag, section_body(id, record));
    // optional sections, ids strictly increasing: the understood panel (I3-D) and any kept verbatim
    std::vector<std::pair<std::uint32_t, const std::string*>> optional;
    std::string panel_body;
    if (record.panel) {
        panel_body = section_body(kPanelMarkerTable, record);
        optional.emplace_back(kPanelMarkerTable, &panel_body);
    }
    for (const auto& [id, body] : record.unknown_optional_sections) optional.emplace_back(id, &body);
    std::sort(optional.begin(), optional.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    for (std::size_t i = 0; i < optional.size(); ++i) {
        if (optional[i].first < kFirstOptionalId || (i > 0 && optional[i].first == optional[i - 1].first)) {
            fail(Kind::invalid_section, "optional section ids must be >= 16 and unique");
        }
        append(optional[i].first, 0, *optional[i].second);
    }
    std::memcpy(&file[0], kMagic, 8);
    put_u32(file, 8, kPreparedLocusSchemaVersion);
    put_u32(file, 12, kHeaderBytes);
    put_u64(file, 16, file.size() - kHeaderBytes);
    put_u32(file, 56, sections);
    put_u32(file, 60, 0);
    reseal_prepared_locus(file);
    return file;
}

void reseal_prepared_locus(std::string& file) {
    if (file.size() < kHeaderBytes) fail(Kind::truncated, "shorter than the 64-byte header");
    const std::array<unsigned char, 32> sum = hex_to_bytes(file_checksum_hex(file));
    std::memcpy(&file[kChecksumAt], sum.data(), kChecksumBytes);
}

namespace {
std::atomic<std::uint64_t>& temp_counter() {
    static std::atomic<std::uint64_t> count{0};
    return count;
}
std::string directory_of(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? std::string(".") : slash == 0 ? std::string("/") : path.substr(0, slash);
}
void write_all(int fd, const std::string& bytes, std::size_t limit, const std::string& tmp) {
    std::size_t done = 0;
    while (done < limit) {
        const ssize_t n = ::write(fd, bytes.data() + done, limit - done);
        if (n < 0) {
            if (errno == EINTR) continue;
            fail(Kind::io, "write failed for " + tmp + ": " + std::strerror(errno));
        }
        done += static_cast<std::size_t>(n);
    }
}
}  // namespace

// Test hook (declared in tests only): 0 = complete publication; 1 = stop after writing half the temporary (an
// interrupted write: the temporary is left, the target never created); 2 = stop after the link, before the unlink;
// I3-D: 3 / 4 / 5 = inject a failure of the temporary's unlink / the directory fsync / the directory open.
void publish_prepared_locus_stage(const std::string& bytes, const std::string& path, int stop_after) {
    const std::string tmp = path + ".tmp." + std::to_string(::getpid()) + "." + std::to_string(++temp_counter());
    const int fd = ::open(tmp.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
    if (fd < 0) fail(Kind::io, "cannot create " + tmp + ": " + std::strerror(errno));
    try {
        if (stop_after == 1) {
            write_all(fd, bytes, bytes.size() / 2, tmp);
            ::close(fd);
            return;
        }
        write_all(fd, bytes, bytes.size(), tmp);
        if (::fsync(fd) != 0) fail(Kind::io, "fsync failed for " + tmp + ": " + std::strerror(errno));
    } catch (...) {
        ::close(fd);
        ::unlink(tmp.c_str());
        throw;
    }
    if (::close(fd) != 0) {
        ::unlink(tmp.c_str());
        fail(Kind::io, "close failed for " + tmp);
    }
    if (::link(tmp.c_str(), path.c_str()) != 0) {
        const int error = errno;
        ::unlink(tmp.c_str());
        if (error == EEXIST) fail(Kind::exists, path + " exists; an index is never overwritten");
        fail(Kind::io, "cannot publish " + path + ": " + std::strerror(error));
    }
    if (stop_after == 2) return;
    // I3-D: after the link the target exists and is complete (its bytes were fsynced before the link). A later failure
    // is reported as post_publication; the valid target is never deleted to pretend nothing was published. Stages 3-5
    // inject a failure of each later step (tests).
    std::string problems;
    errno = 0;
    const int unlinked = stop_after == 3 ? (errno = EACCES, -1) : ::unlink(tmp.c_str());
    if (unlinked != 0) problems += "removing the temporary " + tmp + " failed (" + std::strerror(errno) + "); ";
    const int dir = stop_after == 5 ? (errno = ENOENT, -1) : ::open(directory_of(path).c_str(), O_RDONLY);
    if (dir < 0) {
        problems += "opening the directory to fsync it failed (" + std::string(std::strerror(errno)) + "); ";
    } else {
        const int synced = stop_after == 4 ? (errno = EIO, -1) : ::fsync(dir);
        if (synced != 0) problems += "fsyncing the directory failed (" + std::string(std::strerror(errno)) + "); ";
        ::close(dir);
    }
    if (!problems.empty()) {
        fail(Kind::post_publication, path + " was published (complete; linked after fsync) and is left in place; then " +
                                         problems + "a temporary may remain or the directory entry's durability is unconfirmed");
    }
}

void publish_prepared_locus(const std::string& bytes, const std::string& path) { publish_prepared_locus_stage(bytes, path, 0); }

void write_prepared_locus(const PreparedLocusRecord& record, const std::string& path) {
    publish_prepared_locus(encode_prepared_locus(record), path);
}

PreparedLocusRecord decode_prepared_locus(const std::string& file_bytes, const PreparedLocusLoadLimits& limits,
                                          PreparedLocusLoadReport* report) {
    PreparedLocusLoadReport local;
    if (report == nullptr) report = &local;
    Accounting acc(limits.max_accounted_bytes);
    if (file_bytes.size() > limits.max_file_bytes) fail(Kind::limit_exceeded, "file exceeds max_file_bytes");
    acc.charge(alloc(file_bytes.size()), "read buffer");
    DecodeTallies tallies;
    PreparedLocusRecord r = decode_with(file_bytes, limits, acc, tallies, report);
    if (report != nullptr) {
        report->accounted_peak_bytes = acc.peak();
        report->accounting_limit = acc.limit();
        report->peak_phase = acc.peak_phase();
        report->phase_peaks = acc.phase_peaks();
    }
    return r;
}

namespace {
// Phases 1-3 into `acc`: stat, header, read buffer, checksum, decode, invariants; the buffer released at the end.
PreparedLocusRecord read_record(const std::string& path, const PreparedLocusLoadLimits& limits, Accounting& acc,
                                DecodeTallies& tallies, PreparedLocusLoadReport* report) {
    Clock::time_point mark = Clock::now();
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0) fail(Kind::io, "cannot stat " + path + ": " + std::strerror(errno));
    const std::uint64_t size = static_cast<std::uint64_t>(st.st_size);
    if (size < kHeaderBytes) fail(Kind::truncated, path + " is shorter than the 64-byte header");
    if (size > limits.max_file_bytes) fail(Kind::limit_exceeded, path + " exceeds max_file_bytes");
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) fail(Kind::io, "cannot open " + path + ": " + std::strerror(errno));
    struct Close {
        int fd;
        ~Close() { ::close(fd); }
    } close_fd{fd};
    auto read_exact = [&](char* into, std::size_t n) {
        std::size_t done = 0;
        while (done < n) {
            const ssize_t got = ::read(fd, into + done, n - done);
            if (got < 0 && errno == EINTR) continue;
            if (got <= 0) fail(Kind::truncated, path + " ended early");
            done += static_cast<std::size_t>(got);
        }
    };
    // the header first: nothing is allocated from the payload length until it is shown consistent with the file size
    std::string header(kHeaderBytes, '\0');
    read_exact(&header[0], kHeaderBytes);
    if (std::memcmp(header.data(), kMagic, 8) != 0) fail(Kind::bad_magic, path + " is not a prepared-locus index");
    if (get_u32(header, 8) != kPreparedLocusSchemaVersion) fail(Kind::unsupported_schema, "schema version " + std::to_string(get_u32(header, 8)));
    if (get_u32(header, 12) != kHeaderBytes) fail(Kind::invalid_header, "header_bytes is not 64");
    if (get_u64(header, 16) != size - kHeaderBytes) fail(get_u64(header, 16) > size - kHeaderBytes ? Kind::truncated : Kind::invalid_header,
                                                         "payload_bytes disagrees with the file size");
    acc.charge(alloc(size), "read buffer");
    std::string file;
    file.resize(static_cast<std::size_t>(size));
    std::memcpy(&file[0], header.data(), kHeaderBytes);
    read_exact(&file[kHeaderBytes], static_cast<std::size_t>(size - kHeaderBytes));
    if (report != nullptr) report->seconds.emplace_back("load_read", since(mark));
    PreparedLocusRecord r = decode_with(file, limits, acc, tallies, report);
    std::string().swap(file);
    acc.release(alloc(size));
    return r;
}
}  // namespace

PreparedLocusRecord read_prepared_locus_record(const std::string& path, const PreparedLocusLoadLimits& limits,
                                               PreparedLocusLoadReport* report) {
    PreparedLocusLoadReport local;
    if (report == nullptr) report = &local;
    Accounting acc(limits.max_accounted_bytes);
    DecodeTallies tallies;
    PreparedLocusRecord r = read_record(path, limits, acc, tallies, report);
    if (report != nullptr) {
        report->accounted_peak_bytes = acc.peak();
        report->accounted_final_bytes = acc.used();
        report->accounting_limit = acc.limit();
        report->peak_phase = acc.peak_phase();
        report->phase_peaks = acc.phase_peaks();
    }
    return r;
}

std::shared_ptr<const PreparedLocus> load_prepared_locus(const std::string& path, const PreparedLocusLoadLimits& limits,
                                                         PreparedLocusLoadReport* report) {
    PreparedLocusLoadReport local;
    if (report == nullptr) report = &local;
    Accounting acc(limits.max_accounted_bytes);
    DecodeTallies tallies;
    PreparedLocusRecord r = read_record(path, limits, acc, tallies, report);
    auto timed = [&](const char* name, Clock::time_point mark) {
        if (report != nullptr) report->seconds.emplace_back(name, since(mark));
    };
    // phase 4: the model, projected before the constructor runs
    Clock::time_point mark = Clock::now();
    const ModelProjection projection = project_model(r.input);
    acc.charge(projection.persistent, "model reconstruction (persistent)");
    acc.charge(projection.transient, "model reconstruction (largest block transient)");
    PreparedLocus::Parts parts;
    try {
        parts.model.emplace(prepare_genotype_model(r.input));
    } catch (const std::invalid_argument& e) {
        fail(Kind::invalid_structure, std::string("the model constructor refused the input: ") + e.what());
    }
    acc.release(projection.transient);
    timed("load_model", mark);
    mark = Clock::now();
    parts.fingerprint = model_fingerprint(*parts.model);
    if (parts.fingerprint != r.model_fingerprint) fail(Kind::reconstruction_mismatch, "the model fingerprint differs");
    timed("load_fingerprint", mark);
    mark = Clock::now();
    if (prepared_locus_model_content_sha256(*parts.model) != r.model_content_sha256) {
        fail(Kind::reconstruction_mismatch, "the model-content digest differs");
    }
    timed("load_model_digest", mark);
    // phase 5: the structural source is released (DESIGN B); its decoded charge is exactly the section's tally
    const std::uint64_t input_charge = tallies.input;
    r.input = StructuralLocusInput{};
    acc.release(input_charge);
    // node sequences into the runtime map (strings moved; nodes and buckets charged, the vector released)
    using NodeMapValue = std::pair<const std::string, std::string>;
    acc.charge(add(mul(r.nodes.size(), alloc(16 + sizeof(NodeMapValue))), alloc(mul(8, add(mul(2, r.nodes.size()), 1)))),
               "node map");
    parts.node_sequences.reserve(r.nodes.size());
    for (auto& [id, seq] : r.nodes) parts.node_sequences.emplace(std::move(id), std::move(seq));
    std::vector<std::pair<std::string, std::string>>().swap(r.nodes);
    acc.release(alloc(mul(parts.node_sequences.size(), sizeof(std::pair<std::string, std::string>))));
    PanelBound bound;
    if (r.panel) {
        // I3-D: the persisted panel, already structurally validated against the base (check_panel), verified against
        // the stored complete digest BEFORE any use, then moved into the locus: no pre-count, no build_marker_panel, no
        // second representation. Its storage was charged at decode and stays charged as the retained panel.
        mark = Clock::now();
        if (prepared_locus_panel_content_sha256(*r.panel) != r.panel_content_sha256) {
            fail(Kind::reconstruction_mismatch, "the persisted panel's content digest differs from DERIVED_DIGESTS");
        }
        timed("load_panel_digest", mark);
        parts.panel.emplace(std::move(*r.panel));
        r.panel.reset();
        parts.panel_source = "loaded";
        report->panel_source = "loaded";
        report->panel_decoded_accounted_bytes = tallies.panel;
    } else {
        // the v1 route (DESIGN D phases 6-7): exact pre-count, the panel bound charged, the existing constructor,
        // released to its actual storage, verified against the stored digest
        std::string reason = "no panel section";
        for (const auto& [id, body] : r.unknown_optional_sections) {
            if (id != kPanelMarkerTable) continue;
            reason = body.size() >= 8 ? "unsupported panel section (panel_extension_version " + std::to_string(get_u32(body, 0)) +
                                            ", marker_panel_algorithm_version " + std::to_string(get_u32(body, 4)) + ")"
                                      : "unsupported panel section (shorter than its version fields)";
        }
        mark = Clock::now();
        bound = bound_panel(*parts.model, r.key.kmer_size, r.key.syncmer_s, acc);
        timed("load_panel_precount", mark);
        acc.charge(bound.bytes, "marker panel reconstruction (bound)");
        mark = Clock::now();
        parts.panel.emplace(build_marker_panel(*parts.model, r.key.kmer_size, r.key.syncmer_s));
        acc.release(bound.bytes);
        acc.charge(panel_storage(*parts.panel), "marker panel (retained)");
        timed("load_panel", mark);
        mark = Clock::now();
        if (prepared_locus_panel_content_sha256(*parts.panel) != r.panel_content_sha256) {
            fail(Kind::reconstruction_mismatch, "the panel-content digest differs");
        }
        timed("load_panel_digest", mark);
        parts.panel_source = "rebuilt";
        report->panel_source = "rebuilt";
        report->panel_rebuild_reason = reason;
    }
    parts.key = r.key;
    parts.resolved = r.resolved_reference;
    parts.report = std::move(r.report);
    parts.origin = "loaded";
    parts.file_checksum = report->file_checksum;
    parts.stored_provenance = describe_provenance(r.provenance);
    parts.id = next_prepared_locus_id();
    report->accounted_peak_bytes = acc.peak();
    report->accounted_final_bytes = acc.used();
    report->accounting_limit = acc.limit();
    report->peak_phase = acc.peak_phase();
    report->phase_peaks = acc.phase_peaks();
    report->model_projection_bytes = add(projection.persistent, projection.transient);
    report->panel_projection_bytes = bound.bytes;
    report->panel_syncmer_occurrences = bound.occurrences;
    report->structural_input_accounted_bytes = input_charge;
    ++load_counter();
    return std::make_shared<const PreparedLocus>(std::move(parts));
}

}  // namespace panvar
