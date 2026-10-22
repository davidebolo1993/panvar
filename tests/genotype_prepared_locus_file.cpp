// I3-C: the persistent prepared-locus index v1, focused (no samples are analysed here; that is the e2e gate's job).
//
//   genotype_prepared_locus_file_test SCRATCH_DIR
//
// A locus built by the REAL builder from a hand-written graph with a curated partial template, mixed template
// orientation, alias walks (two walks, one sequence) and an empty allele; then: independent preparations give identical
// bytes; write -> read -> rewrite is byte-identical and keeps the stored provenance; load reconstructs the same model,
// panel, fingerprint and report; versions and provenance; header and payload corruption (truncation, bit flips, semantic
// damage behind a recomputed checksum); reconstruction mismatch; cumulative allocation accounting controls at each phase;
// publication refusal and interrupted writes.
#include "panvar/genotype_command.hpp"
#include "panvar/prepared_locus.hpp"
#include "panvar/prepared_locus_file.hpp"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
// test hook (prepared_locus_file.cpp): 1 = stop after half the temporary is written; 2 = stop after the link
void publish_prepared_locus_stage(const std::string& bytes, const std::string& path, int stop_after);
}
#include "panvar/marker_model.hpp"
#include "panvar/syncmer.hpp"

namespace {

using Kind = panvar::PreparedLocusFormatError::Kind;
std::string g_dir;
std::string path(const std::string& name) { return g_dir + "/" + name; }

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string random_bases(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}

std::string slurp(const std::string& file) {
    std::ifstream in(file, std::ios::binary);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

void remove_exact(const std::string& file) {   // never a wildcard: one named file
    std::error_code ignored;
    std::filesystem::remove(file, ignored);
}

// Expects `f` to refuse with `kind`; returns the message.
std::string expect(Kind kind, const std::function<void()>& f, const std::string& what) {
    try {
        f();
    } catch (const panvar::PreparedLocusFormatError& e) {
        require(e.kind() == kind, what + ": expected " + panvar::prepared_locus_error_name(kind) + ", got: " + e.what());
        return e.what();
    }
    throw std::runtime_error(what + ": expected " + std::string(panvar::prepared_locus_error_name(kind)) + ", accepted");
}

// The hand-written graph. Bubbles: SV1 (node 2 vs 3 vs 3b, where 3b spells 3's sequence: an alias walk), SV2 (node 5
// vs the direct 4->6 edge: an empty allele), SV3 (7 vs 8). Templates: A (reference), B, C (the alias), D written in
// reverse (a Reverse-oriented template), E stops after SV2 (coverage-partial; curated out by the manifest).
void write_fixture_graph(const std::string& file, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    const std::string n1 = random_bases(rng, 400), n2 = random_bases(rng, 220), n3 = random_bases(rng, 70), n4 = random_bases(rng, 260),
                      n5 = random_bases(rng, 160), n6 = random_bases(rng, 260), n7 = random_bases(rng, 120),
                      n8 = random_bases(rng, 40), n9 = random_bases(rng, 400);
    std::ofstream out(file);
    out << "H\tVN:Z:1.0\n";
    const std::vector<std::pair<std::string, std::string>> nodes = {{"1", n1}, {"2", n2}, {"3", n3}, {"3b", n3}, {"4", n4},
                                                                    {"5", n5}, {"6", n6}, {"7", n7}, {"8", n8}, {"9", n9}};
    for (const auto& [id, seq] : nodes) out << "S\t" << id << '\t' << seq << '\n';
    for (const char* edge : {"1\t+\t2\t+", "1\t+\t3\t+", "1\t+\t3b\t+", "2\t+\t4\t+", "3\t+\t4\t+", "3b\t+\t4\t+", "4\t+\t5\t+",
                             "5\t+\t6\t+", "4\t+\t6\t+", "6\t+\t7\t+", "6\t+\t8\t+", "7\t+\t9\t+", "8\t+\t9\t+"}) {
        out << "L\t" << edge << "\t0M\n";
    }
    out << "P\tsampleA#1#chr1:0-100\t1+,2+,4+,5+,6+,7+,9+\t*\n";
    out << "P\tsampleB#1#chr1:0-100\t1+,3+,4+,6+,8+,9+\t*\n";
    out << "P\tsampleC#1#chr1:0-100\t1+,3b+,4+,5+,6+,8+,9+\t*\n";
    out << "P\tsampleD#1#chr1:0-100\t9-,7-,6-,4-,2-,1-\t*\n";
    out << "P\tsampleE#1#chr1:0-100\t1+,2+,4+,5+,6+\t*\n";
}

panvar::GenotypeCommandOptions fixture_options(const std::string& gfa, const std::string& prefix) {
    panvar::GenotypeCommandOptions o;
    o.gfa_path = path(gfa);
    o.reference_path = "sampleA#1#chr1:0-100";
    o.curation_manifest_path = path("fixture_curation.txt");
    o.out_prefix = path(prefix);
    o.kmer_size = 21;
    o.syncmer_s = 11;
    o.min_variant_bp = 30;
    return o;
}

// The body offset and length of section `id` in a file.
std::pair<std::size_t, std::size_t> section(const std::string& file, std::uint32_t id) {
    std::size_t at = 64;
    while (at + 16 <= file.size()) {
        std::uint32_t sid = 0;
        std::uint64_t bytes = 0;
        std::memcpy(&sid, &file[at], 4);
        std::memcpy(&bytes, &file[at + 8], 8);
        if (sid == id) return {at + 16, static_cast<std::size_t>(bytes)};
        at += 16 + static_cast<std::size_t>(bytes);
    }
    throw std::runtime_error("section not found");
}
void put_u64(std::string& s, std::size_t at, std::uint64_t v) { std::memcpy(&s[at], &v, 8); }
void put_u32(std::string& s, std::size_t at, std::uint32_t v) { std::memcpy(&s[at], &v, 4); }
std::uint64_t get_u64(const std::string& s, std::size_t at) {
    std::uint64_t v = 0;
    std::memcpy(&v, &s[at], 8);
    return v;
}

const panvar::PreparedLocusLoadLimits kDefault{};

panvar::PreparedLocusRecord decode(const std::string& bytes, const panvar::PreparedLocusLoadLimits& limits = kDefault) {
    return panvar::decode_prepared_locus(bytes, limits);
}

std::shared_ptr<const panvar::PreparedLocus> load_bytes(const std::string& bytes, const std::string& name,
                                                        const panvar::PreparedLocusLoadLimits& limits = kDefault,
                                                        panvar::PreparedLocusLoadReport* report = nullptr) {
    remove_exact(path(name));
    panvar::publish_prepared_locus(bytes, path(name));
    auto locus = panvar::load_prepared_locus(path(name), limits, report);
    remove_exact(path(name));
    return locus;
}

// ---------------------------------------------------------------------------------------------------------------------
struct Fixture {
    panvar::PreparedLocusPreparation prepared;
    std::string bytes;
};

Fixture build_fixture() {
    write_fixture_graph(path("fixture.gfa"), 11);
    {
        std::ofstream manifest(path("fixture_curation.txt"));
        manifest << "# the partial template\nsampleE#1#chr1:0-100\n";
    }
    Fixture f;
    f.prepared = panvar::prepare_genotype_locus_record(fixture_options("fixture.gfa", "fixture_prep1"));
    f.bytes = panvar::encode_prepared_locus(f.prepared.record);
    const panvar::PreparedLocusRecord& r = f.prepared.record;
    // the fixture really exercises the builder semantics the loader must preserve
    bool curated_partial = false, curated_in_coverage = false, alias = false, empty = false;
    for (const auto& e : r.report.curated_exclusions) {
        curated_partial |= e.template_name == "sampleE#1#chr1:0-100" && e.reason.find("locus_coverage_partial") != std::string::npos;
    }
    for (const auto& c : r.report.non_traversing_coverage) curated_in_coverage |= c.template_name == "sampleE#1#chr1:0-100";
    for (const auto& b : r.report.blocks) alias |= b.canonical_alleles < b.source_alleles;
    for (const auto& b : r.input.blocks) {
        for (const auto& a : b.alleles) empty |= a.source_id == "empty";
    }
    std::size_t forward = 0, reverse = 0;
    for (const auto& t : r.input.templates) {
        (t.source_orientation == panvar::ChainOrientation::Forward ? forward : reverse) += 1;
        require(t.name != "sampleE#1#chr1:0-100", "fixture: the curated template must not be an input template");
    }
    require(curated_partial, "fixture: the partial template must be curated with its partial reason");
    require(curated_in_coverage, "fixture: the curated template must still appear in the coverage table");
    require(alias, "fixture: an alias walk must collapse into one canonical allele");
    require(empty, "fixture: the builder's empty-interior source must be present");
    require(forward > 0 && reverse > 0, "fixture: template orientations must be mixed");
    require(r.report.forward_templates + r.report.reverse_templates == r.input.templates.size() + 1,
            "fixture: orientation counts include the curated template (counted before curation)");
    require(r.report.refused_templates.empty() && r.report.partial_templates.empty(), "fixture: curation removes the refusal records");
    return f;
}

void test_determinism_round_trip_rewrite(const Fixture& f) {
    // independent preparations: identical bytes
    const panvar::PreparedLocusPreparation again = panvar::prepare_genotype_locus_record(fixture_options("fixture.gfa", "fixture_prep2"));
    require(panvar::encode_prepared_locus(again.record) == f.bytes, "two independent preparations must give identical files");
    require(again.locus->id() != f.prepared.locus->id(), "with different runtime ids: no runtime id is written");
    // prepare -> write -> read -> rewrite
    remove_exact(path("fixture.pvpl"));
    remove_exact(path("fixture_rewrite.pvpl"));
    panvar::write_prepared_locus(f.prepared.record, path("fixture.pvpl"));
    require(slurp(path("fixture.pvpl")) == f.bytes, "the written file is the encoding");
    const panvar::PreparedLocusRecord read = panvar::read_prepared_locus_record(path("fixture.pvpl"), kDefault);
    require(panvar::encode_prepared_locus(read) == f.bytes, "decode then encode is the identity");
    panvar::write_prepared_locus(read, path("fixture_rewrite.pvpl"));
    require(slurp(path("fixture_rewrite.pvpl")) == f.bytes, "rewriting an unchanged index gives identical bytes");
    // load: the reconstructed object equals the prepared one
    panvar::PreparedLocusLoadReport report;
    const std::uint64_t loads = panvar::prepared_locus_loads(), builds = panvar::prepared_locus_builds();
    const auto loaded = panvar::load_prepared_locus(path("fixture.pvpl"), kDefault, &report);
    require(panvar::prepared_locus_loads() == loads + 1 && panvar::prepared_locus_builds() == builds, "a load is not a preparation");
    const panvar::PreparedLocus& p = *f.prepared.locus;
    require(loaded->fingerprint() == p.fingerprint(), "fingerprint");
    require(panvar::prepared_locus_model_content_sha256(loaded->model()) == panvar::prepared_locus_model_content_sha256(p.model()),
            "model content");
    require(panvar::prepared_locus_panel_content_sha256(loaded->panel()) == panvar::prepared_locus_panel_content_sha256(p.panel()),
            "panel content");
    require(loaded->node_sequences() == p.node_sequences() && loaded->resolved_reference() == p.resolved_reference(), "nodes, reference");
    require(loaded->key().describe() == p.key().describe(), "identity key");
    panvar::PreparedLocusRecord with_loaded_report = f.prepared.record;
    with_loaded_report.report = loaded->report();
    require(panvar::encode_prepared_locus(with_loaded_report) == f.bytes, "the loaded report equals the prepared report");
    require(loaded->digest() == p.digest(), "the I3-A convenience digest agrees too");
    require(loaded->origin() == "loaded" && p.origin() == "prepared" && loaded->file_checksum() == report.file_checksum &&
                loaded->id() != p.id(),
            "runtime facts: origin, file checksum, a fresh runtime id");
    require(report.structural_input_accounted_bytes > 0 && report.accounted_final_bytes < report.accounted_peak_bytes,
            "the structural input is released after reconstruction");
    remove_exact(path("fixture.pvpl"));
    remove_exact(path("fixture_rewrite.pvpl"));
    std::cout << "  ok: independent preparations identical; write/read/rewrite byte-identical; load reconstructs model, panel, "
                 "fingerprint, report, nodes (" << f.bytes.size() << " bytes)\n";
}

// An allele whose sequence is empty (the builder's divergence trim keeps boundary bases, so the fixture graph cannot give
// one; the model constructor accepts it): added to a variant block with the report kept consistent and the derived
// digests recomputed through the public constructors. It must round-trip byte-identically and load.
void test_empty_sequence_allele(const Fixture& f) {
    panvar::PreparedLocusRecord r = f.prepared.record;
    std::size_t bi = 0;
    while (r.input.blocks[bi].kind != panvar::GenotypeBlockKind::Variant) ++bi;
    r.input.blocks[bi].alleles.push_back({"bypass", ""});
    panvar::LocusBlockReport& br = r.report.blocks[bi];
    br.source_alleles += 1;
    br.canonical_alleles += 1;
    br.min_bp = 0;
    r.report.total_source_alleles += 1;
    r.report.total_canonical_alleles += 1;
    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(r.input);
    r.model_fingerprint = panvar::model_fingerprint(model);
    r.model_content_sha256 = panvar::prepared_locus_model_content_sha256(model);
    r.panel_content_sha256 = panvar::prepared_locus_panel_content_sha256(panvar::build_marker_panel(model, 21, 11));
    const std::string bytes = panvar::encode_prepared_locus(r);
    require(panvar::encode_prepared_locus(decode(bytes)) == bytes, "an empty-sequence allele round-trips");
    const auto loaded = load_bytes(bytes, "empty_allele.pvpl");
    bool found = false;
    for (const auto& a : loaded->model().blocks()[bi].alleles) found |= a.sequence.empty();
    require(found && loaded->fingerprint() == r.model_fingerprint, "and loads with its empty allele");
    std::cout << "  ok: an empty-sequence allele round-trips and loads\n";
}

void test_versions_and_provenance(const Fixture& f) {
    // a different (compatible) build: accepted; its provenance is reported and preserved by a rewrite
    panvar::PreparedLocusRecord other = f.prepared.record;
    other.provenance.compiler = "some other compiler 0.0";
    other.provenance.source_commit = "0000000000000000000000000000000000000000";
    const std::string other_bytes = panvar::encode_prepared_locus(other);
    const auto loaded = load_bytes(other_bytes, "prov.pvpl");
    require(loaded->stored_provenance().find("some other compiler 0.0") != std::string::npos, "the stored provenance is reported");
    require(loaded->fingerprint() == f.prepared.locus->fingerprint(), "provenance never gates loading");
    require(panvar::encode_prepared_locus(decode(other_bytes)) == other_bytes, "a rewrite keeps the stored provenance, not this build's");
    require(!(panvar::current_prepared_locus_provenance() == other.provenance), "fixture: the provenances differ");
    // an incompatible preparation algorithm
    panvar::PreparedLocusRecord future = f.prepared.record;
    future.preparation_algorithm_version = panvar::kPreparationAlgorithmVersion + 1;
    expect(Kind::incompatible_preparation, [&] { decode(panvar::encode_prepared_locus(future)); }, "preparation version");
    // an unsupported schema
    std::string schema = f.bytes;
    put_u32(schema, 8, 2);
    panvar::reseal_prepared_locus(schema);
    expect(Kind::unsupported_schema, [&] { decode(schema); }, "schema version");
    std::cout << "  ok: a different compatible build accepted with its provenance reported and preserved by rewrite; "
                 "preparation version and schema refused\n";
}

void test_container_corruption(const Fixture& f) {
    const std::string& b = f.bytes;
    // truncation at every section boundary, inside the header, and at deterministic pseudo-random offsets
    std::vector<std::size_t> cuts = {0, 1, 8, 63, 64, 65};
    for (std::uint32_t id = 1; id <= 7; ++id) {
        const auto [at, len] = section(b, id);
        cuts.push_back(at - 16);
        cuts.push_back(at);
        cuts.push_back(at + len / 2);
    }
    std::mt19937_64 rng(5);
    for (int i = 0; i < 40; ++i) cuts.push_back(rng() % b.size());
    for (const std::size_t cut : cuts) {
        if (cut >= b.size()) continue;
        const std::string t = b.substr(0, cut);
        try {
            decode(t);
            throw std::runtime_error("a truncation was accepted at " + std::to_string(cut));
        } catch (const panvar::PreparedLocusFormatError& e) {
            require(e.kind() == Kind::truncated || e.kind() == Kind::bad_magic,
                    "truncation at " + std::to_string(cut) + ": " + e.what());
        }
    }
    // a truncated FILE is refused by the bounded reader before anything is allocated from its claimed payload
    remove_exact(path("trunc.pvpl"));
    {
        std::ofstream out(path("trunc.pvpl"), std::ios::binary);
        out << b.substr(0, b.size() - 7);
    }
    expect(Kind::truncated, [&] { panvar::read_prepared_locus_record(path("trunc.pvpl"), kDefault); }, "truncated file");
    remove_exact(path("trunc.pvpl"));
    // single-bit flips: header fields refuse with their own kind; the checksum covers every other header byte and the
    // whole payload (section_count and flags included)
    for (std::size_t at = 0; at < 64; ++at) {
        std::string x = b;
        x[at] = static_cast<char>(x[at] ^ 0x10);
        try {
            decode(x);
            throw std::runtime_error("a header bit flip was accepted at " + std::to_string(at));
        } catch (const panvar::PreparedLocusFormatError& e) {
            const Kind expected = at < 8 ? Kind::bad_magic : at < 12 ? Kind::unsupported_schema : at < 16 ? Kind::invalid_header
                                : at < 24 ? Kind::truncated : Kind::checksum_mismatch;
            require(e.kind() == expected || (at >= 16 && at < 24 && e.kind() == Kind::invalid_header),
                    "header flip at " + std::to_string(at) + ": " + e.what());
        }
    }
    for (std::uint32_t id = 1; id <= 7; ++id) {
        const auto [at, len] = section(b, id);
        for (const std::size_t off : {at - 16, at - 8, at, at + len / 3, at + len - 1}) {
            std::string x = b;
            x[off] = static_cast<char>(x[off] ^ 0x01);
            expect(Kind::checksum_mismatch, [&] { decode(x); }, "payload flip in section " + std::to_string(id));
        }
    }
    // semantic header damage behind a recomputed checksum
    {
        std::string x = b;
        put_u32(x, 60, 4);
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_header, [&] { decode(x); }, "reserved header flags");
        std::string y = b;
        put_u32(y, 56, 8);
        panvar::reseal_prepared_locus(y);
        expect(Kind::invalid_section, [&] { decode(y); }, "section_count");
    }
    std::cout << "  ok: " << cuts.size() << " truncations; every header byte; 35 payload bit flips (checksum covers "
                 "section_count and flags); resealed flags and section_count refused\n";
}

// Primitive-level damage at known offsets, resealed.
void test_semantic_bytes(const Fixture& f) {
    const std::string& b = f.bytes;
    const auto [id_at, id_len] = section(b, 1);
    (void)id_len;
    {   // a string length beyond the section
        std::string x = b;
        put_u64(x, id_at, 1ull << 40);
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_length, [&] { decode(x); }, "string length past the end");
    }
    const auto [in_at, in_len] = section(b, 4);
    (void)in_len;
    const std::uint64_t locus_len = get_u64(b, in_at);
    const std::size_t orientation_at = in_at + 8 + static_cast<std::size_t>(locus_len), blocks_at = orientation_at + 1;
    {   // an invalid enum
        std::string x = b;
        x[orientation_at] = 7;
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_enum, [&] { decode(x); }, "enum");
    }
    {   // a count of 2^63: refused by the element limit; with the limit lifted, by the remaining-bytes check
        std::string x = b;
        put_u64(x, blocks_at, 1ull << 63);
        panvar::reseal_prepared_locus(x);
        expect(Kind::limit_exceeded, [&] { decode(x); }, "count 2^63, element limit");
        panvar::PreparedLocusLoadLimits open = kDefault;
        open.max_elements = ~0ull;
        expect(Kind::overflow, [&] { decode(x, open); }, "count 2^63 x element size overflows");
        put_u64(x, blocks_at, 1ull << 20);
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_count, [&] { decode(x, open); }, "a count that cannot fit");
    }
    const auto [co_at, co_len] = section(b, 2);
    (void)co_len;
    {   // an invalid bool
        std::string x = b;
        x[co_at + 4] = 2;
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_enum, [&] { decode(x); }, "bool");
    }
    {   // a section length that does not match its content: one trailing byte inside COMPATIBILITY
        std::string x = b;
        x.insert(co_at + 6, 1, '\0');
        put_u64(x, co_at - 8, get_u64(x, co_at - 8) + 1);
        put_u64(x, 16, get_u64(x, 16) + 1);
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_length, [&] { decode(x); }, "section length longer than its content");
        std::string y = b;
        y.erase(co_at + 5, 1);
        put_u64(y, co_at - 8, get_u64(y, co_at - 8) - 1);
        put_u64(y, 16, get_u64(y, 16) - 1);
        panvar::reseal_prepared_locus(y);
        expect(Kind::invalid_length, [&] { decode(y); }, "section length shorter than its content");
    }
    // sections: unknown optional (skipped), unknown required, reserved id, reserved flag, order, missing
    auto with_section = [&](std::uint32_t id, std::uint32_t flags, const std::string& body) {
        std::string x = b;
        std::string head(16, '\0');
        put_u32(head, 0, id);
        put_u32(head, 4, flags);
        put_u64(head, 8, body.size());
        x += head + body;
        put_u64(x, 16, x.size() - 64);
        put_u32(x, 56, 8);
        panvar::reseal_prepared_locus(x);
        return x;
    };
    {
        const std::string optional = with_section(40, 0, "future extension");
        const auto loaded = load_bytes(optional, "optional.pvpl");
        require(loaded->fingerprint() == f.prepared.locus->fingerprint(), "an unknown optional section is skipped");
        expect(Kind::unknown_required_section, [&] { decode(with_section(40, 1, "x")); }, "unknown required section");
        expect(Kind::invalid_section, [&] { decode(with_section(9, 0, "x")); }, "reserved section id");
        expect(Kind::invalid_section, [&] { decode(with_section(40, 2, "x")); }, "reserved section flag");
        expect(Kind::invalid_section, [&] { decode(with_section(3, 1, "x")); }, "a repeated (non-increasing) section id");
        std::string x = b;   // a required section flagged optional
        put_u32(x, co_at - 12, 0);
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_section, [&] { decode(x); }, "a required section flagged optional");
    }
    std::cout << "  ok: string length, enum, bool, count limit / overflow / fit, section length +-1, unknown optional "
                 "skipped, unknown required, reserved id and flag, repeated id, required-flag refusals\n";
}

// Semantic damage at the record level: the encoder recomputes the checksum, the invariants or digests must refuse.
void test_semantic_records(const Fixture& f) {
    const panvar::PreparedLocusRecord& base = f.prepared.record;
    struct Case {
        std::string what;
        std::function<void(panvar::PreparedLocusRecord&)> change;
    };
    const std::vector<Case> invariants = {
        {"duplicate block id", [](auto& r) { r.input.blocks[1].block_id = r.input.blocks[0].block_id; }},
        {"duplicate source", [](auto& r) {
             for (auto& b : r.input.blocks) if (b.alleles.size() > 1) { b.alleles[1].source_id = b.alleles[0].source_id; return; }
         }},
        {"duplicate template", [](auto& r) { r.input.templates[1].name = r.input.templates[0].name; }},
        {"duplicate node", [](auto& r) { r.nodes[1].first = r.nodes[0].first; }},
        {"unsorted nodes", [](auto& r) { std::swap(r.nodes[0], r.nodes[1]); }},
        {"template names a missing source", [](auto& r) { r.input.templates[0].allele_source_ids[0] = "walk:nowhere+"; }},
        {"walk names a missing node", [](auto& r) {
             for (auto& b : r.input.blocks) for (auto& a : b.alleles) if (a.source_id.rfind("walk:", 0) == 0) {
                 for (auto& t : r.input.templates) for (auto& s : t.allele_source_ids) if (s == a.source_id) s = "walk:999+";
                 a.source_id = "walk:999+";
                 return;
             }
         }},
        {"malformed walk", [](auto& r) {
             for (auto& b : r.input.blocks) for (auto& a : b.alleles) if (a.source_id.rfind("walk:", 0) == 0) {
                 for (auto& t : r.input.templates) for (auto& s : t.allele_source_ids) if (s == a.source_id) s = "walk:4";
                 a.source_id = "walk:4";
                 return;
             }
         }},
        {"wrong invariant segment count", [](auto& r) { r.input.invariant_segments.push_back(""); }},
        {"report block counts", [](auto& r) { r.report.variant_blocks += 1; }},
        {"report allele totals", [](auto& r) { r.report.total_source_alleles += 1; }},
        {"block report canonical count", [](auto& r) { r.report.blocks[0].canonical_alleles += 1; }},
        {"block report length range", [](auto& r) { r.report.blocks[0].max_bp += 1; }},
        {"block report kind", [](auto& r) {
             r.report.blocks[0].kind = r.report.blocks[0].kind == panvar::GenotypeBlockKind::Variant ? panvar::GenotypeBlockKind::Context
                                                                                                      : panvar::GenotypeBlockKind::Variant;
         }},
        {"orientation counts below the complete templates", [](auto& r) { r.report.forward_templates = 0; r.report.reverse_templates = 0; }},
        {"orientation counts above input + curated", [](auto& r) { r.report.forward_templates += 5; }},
        {"curated name is an input template", [](auto& r) { r.report.curated_exclusions[0].template_name = r.input.templates[0].name; }},
        {"curation digest without exclusions", [](auto& r) { r.report.curated_exclusions.clear(); }},
        {"coverage names an unknown template", [](auto& r) { r.report.non_traversing_coverage[0].template_name = "nobody"; }},
        {"coverage names an unknown block", [](auto& r) { r.report.non_traversing_coverage[0].block_id = "sv:99"; }},
        {"refused without a partial record", [](auto& r) { r.report.refused_templates.push_back(r.input.templates[0].name); }},
        {"representation-dependent block", [](auto& r) {
             r.report.blocks[0].representation_dependent = true;
             r.report.representation_dependent_blocks.push_back(r.report.blocks[0].block_id);
         }},
        {"resolved reference differs", [](auto& r) { r.resolved_reference += "x"; }},
        {"preparation checks not passed", [](auto& r) { r.round_trip_passed = false; }},
        {"boundary node absent", [](auto& r) { r.report.blocks[0].source_node = "nope"; }},
    };
    for (const Case& c : invariants) {
        panvar::PreparedLocusRecord r = base;
        c.change(r);
        expect(Kind::invalid_structure, [&] { decode(panvar::encode_prepared_locus(r)); }, c.what);
    }
    // an uncurated partial template decodes (the builder writes it so) and is refused by the model constructor on load
    {
        panvar::PreparedLocusRecord r = base;
        panvar::PanelTemplateInput partial;
        partial.name = "sampleE#1#chr1:0-100";
        partial.frame_status = panvar::TemplateFrameStatus::Partial;
        partial.partial_reason = "locus_coverage_partial: missing right-terminal coverage";
        r.input.templates.push_back(partial);
        r.report.curated_exclusions.clear();
        r.report.curation_digest.clear();
        r.report.refused_templates.push_back(partial.name);
        panvar::PartialTemplateReport report;
        report.name = partial.name;
        report.reason = partial.partial_reason;
        report.missing_terminal_anchor = "right";
        r.report.partial_templates.push_back(report);
        const std::string bytes = panvar::encode_prepared_locus(r);
        (void)decode(bytes);
        expect(Kind::invalid_structure, [&] { load_bytes(bytes, "partial.pvpl"); }, "an uncurated partial template");
    }
    // digests: altered fingerprint or content digest, or content changed behind unchanged digests
    const std::vector<Case> digests = {
        {"fingerprint altered", [](auto& r) { r.model_fingerprint[0] = r.model_fingerprint[0] == '0' ? '1' : '0'; }},
        {"model digest altered", [](auto& r) { r.model_content_sha256[0] = r.model_content_sha256[0] == '0' ? '1' : '0'; }},
        {"panel digest altered", [](auto& r) { r.panel_content_sha256[0] = r.panel_content_sha256[0] == '0' ? '1' : '0'; }},
        {"an allele sequence changed (lengths kept)", [](auto& r) {
             for (auto& b : r.input.blocks) for (auto& a : b.alleles) if (!a.sequence.empty()) {
                 a.sequence[0] = a.sequence[0] == 'A' ? 'C' : 'A';
                 return;
             }
         }},
    };
    for (const Case& c : digests) {
        panvar::PreparedLocusRecord r = base;
        c.change(r);
        const std::string bytes = panvar::encode_prepared_locus(r);
        (void)decode(bytes);   // structurally valid
        expect(Kind::reconstruction_mismatch, [&] { load_bytes(bytes, "mismatch.pvpl"); }, c.what);
    }
    std::cout << "  ok: " << invariants.size() << " invariant violations refused behind valid checksums; an uncurated partial "
              << "template refused by the model constructor; " << digests.size() << " reconstruction mismatches refused\n";
}

void test_accounting(const Fixture& f) {
    // phase by phase, from a successful load's accounting: a limit one byte below each new running maximum refuses at
    // exactly that phase
    panvar::PreparedLocusLoadReport report;
    (void)load_bytes(f.bytes, "acct.pvpl", kDefault, &report);
    require(report.accounted_peak_bytes > f.bytes.size(), "accounting covers more than the file's bytes");
    require(report.panel_projection_bytes > 0 && report.model_projection_bytes > 0 && report.panel_syncmer_occurrences > 0,
            "the model and panel reconstructions are projected");
    std::uint64_t running = 0;
    std::vector<std::string> refused_at;
    bool model = false, panel = false, buffer = false, decoded = false;
    for (const auto& [label, peak] : report.phase_peaks) {
        if (peak <= running) continue;
        running = peak;
        panvar::PreparedLocusLoadLimits tight = kDefault;
        tight.max_accounted_bytes = peak - 1;
        const std::string message = expect(Kind::limit_exceeded, [&] { load_bytes(f.bytes, "acct.pvpl", tight); }, "limit at " + label);
        require(message.find(label) != std::string::npos, "the refusal names its phase: " + label + " / " + message);
        refused_at.push_back(label);
        buffer |= label == "read buffer";
        model |= label.find("model reconstruction") != std::string::npos;
        panel |= label.find("marker panel reconstruction") != std::string::npos;
        decoded |= label == "allele sequence" || label == "nodes" || label == "blocks";
    }
    require(buffer && panel, "the read buffer and the panel bound are both enforced phases");
    (void)model;
    (void)decoded;
    // container amplification: many tiny nodes whose encoded size passes every byte check but whose decoded containers
    // exceed the accounting limit
    {
        panvar::PreparedLocusRecord r = f.prepared.record;
        for (int i = 0; i < 200000; ++i) {
            char id[16];
            std::snprintf(id, sizeof id, "zz%07d", i);
            r.nodes.emplace_back(id, "");
        }
        const std::string bytes = panvar::encode_prepared_locus(r);
        panvar::PreparedLocusLoadLimits tight = kDefault;
        tight.max_accounted_bytes = 2 * bytes.size();   // more than the file and its read buffer together
        const std::string message = expect(Kind::limit_exceeded, [&] { decode(bytes, tight); }, "node container amplification");
        require(message.find("nodes") != std::string::npos || message.find("node id") != std::string::npos,
                "the amplification is caught while decoding nodes: " + message);
    }
    // nested counts: many template source ids (empty strings, 8 encoded bytes, 24+ decoded)
    {
        panvar::PreparedLocusRecord r = f.prepared.record;
        r.input.templates[0].allele_source_ids.resize(400000);
        const std::string bytes = panvar::encode_prepared_locus(r);
        panvar::PreparedLocusLoadLimits tight = kDefault;
        tight.max_accounted_bytes = 2 * bytes.size();
        expect(Kind::limit_exceeded, [&] { decode(bytes, tight); }, "nested template source-id amplification");
    }
    // the file-size and element limits
    {
        panvar::PreparedLocusLoadLimits tiny = kDefault;
        tiny.max_file_bytes = f.bytes.size() - 1;
        expect(Kind::limit_exceeded, [&] { load_bytes(f.bytes, "acct.pvpl", tiny); }, "max_file_bytes");
        tiny = kDefault;
        tiny.max_elements = 2;
        expect(Kind::limit_exceeded, [&] { decode(f.bytes, tiny); }, "max_elements");
    }
    std::cout << "  ok: accounting peak " << report.accounted_peak_bytes << " B (panel bound " << report.panel_projection_bytes
              << " B, model " << report.model_projection_bytes << " B) for a " << f.bytes.size() << "-byte file; refused one byte "
              << "below each of " << refused_at.size() << " phase maxima; node and nested-count amplification refused; "
                 "file and element limits\n";
}

void test_publication(const Fixture& f) {
    const std::string target = path("pub.pvpl");
    remove_exact(target);
    panvar::publish_prepared_locus(f.bytes, target);
    // an existing target is never overwritten
    panvar::PreparedLocusRecord other = f.prepared.record;
    other.provenance.compiler = "different";
    expect(Kind::exists, [&] { panvar::write_prepared_locus(other, target); }, "existing target");
    require(slurp(target) == f.bytes, "the existing target is untouched");
    // a missing directory: nothing is created
    expect(Kind::io, [&] { panvar::publish_prepared_locus(f.bytes, path("no_such_dir/pub.pvpl")); }, "missing directory");
    require(!std::filesystem::exists(path("no_such_dir")), "nothing is created for a failed publication");
    remove_exact(target);
    // an interrupted write (half the temporary): the target never appears; the loader is unaffected; a retry publishes
    const auto temporaries = [&] {
        std::vector<std::string> out;
        for (const auto& entry : std::filesystem::directory_iterator(g_dir)) {
            const std::string name = entry.path().filename().string();
            if (name.rfind("pub.pvpl.tmp.", 0) == 0) out.push_back(entry.path().string());
        }
        return out;
    };
    for (const std::string& stale : temporaries()) remove_exact(stale);   // earlier runs of this test, by exact name
    panvar::publish_prepared_locus_stage(f.bytes, target, 1);
    require(!std::filesystem::exists(target), "an interrupted write leaves no target");
    require(temporaries().size() == 1 && std::filesystem::file_size(temporaries()[0]) == f.bytes.size() / 2,
            "an interrupted write leaves exactly its partial temporary");
    expect(Kind::io, [&] { panvar::load_prepared_locus(target, kDefault); }, "the loader never opens a temporary");
    panvar::publish_prepared_locus(f.bytes, target);
    require(slurp(target) == f.bytes, "a later publication succeeds beside the stale temporary");
    remove_exact(target);
    for (const std::string& stale : temporaries()) remove_exact(stale);
    // interrupted after the link: the target is complete and loads
    panvar::publish_prepared_locus_stage(f.bytes, target, 2);
    require(slurp(target) == f.bytes && temporaries().size() == 1, "after the link: a complete target and a leftover temporary");
    require(panvar::load_prepared_locus(target, kDefault)->fingerprint() == f.prepared.locus->fingerprint(), "and it loads");
    remove_exact(target);
    for (const std::string& stale : temporaries()) remove_exact(stale);
    std::cout << "  ok: no overwrite (target untouched); missing directory; interrupted write leaves only its temporary and "
                 "no target; interruption after the link leaves a complete loadable target\n";
}


// ---- I3-D: the persisted marker-panel extension ----------------------------------------------------------------------
struct PanelFixture {
    panvar::PreparedLocusPreparation prepared;
    std::string bytes;
};

PanelFixture build_panel_fixture(const Fixture& f) {
    PanelFixture p;
    p.prepared = panvar::prepare_genotype_locus_record(fixture_options("fixture.gfa", "fixture_panel1"), true);
    p.bytes = panvar::encode_prepared_locus(p.prepared.record);
    require(p.prepared.record.panel.has_value(), "the record carries the panel when asked");
    panvar::PreparedLocusRecord without = p.prepared.record;
    without.panel.reset();
    require(panvar::encode_prepared_locus(without) == f.bytes, "without the panel the file is byte-identical to v1");
    return p;
}

void test_panel_round_trip(const Fixture& f, const PanelFixture& p) {
    const auto again = panvar::prepare_genotype_locus_record(fixture_options("fixture.gfa", "fixture_panel2"), true);
    require(panvar::encode_prepared_locus(again.record) == p.bytes, "independent panel preparations write identical bytes");
    require(panvar::encode_prepared_locus(decode(p.bytes)) == p.bytes, "panel extension: decode then encode is the identity");
    remove_exact(path("panel.pvpl"));
    remove_exact(path("panel_rewrite.pvpl"));
    panvar::write_prepared_locus(p.prepared.record, path("panel.pvpl"));
    panvar::write_prepared_locus(panvar::read_prepared_locus_record(path("panel.pvpl"), kDefault), path("panel_rewrite.pvpl"));
    require(slurp(path("panel_rewrite.pvpl")) == p.bytes, "panel extension: rewrite byte-identical");
    // the loaded-panel route builds nothing and pre-counts nothing
    const std::uint64_t builds = panvar::marker_panel_builds(), precounts = panvar::prepared_locus_panel_precounts();
    panvar::PreparedLocusLoadReport report;
    const auto loaded = panvar::load_prepared_locus(path("panel.pvpl"), kDefault, &report);
    require(panvar::marker_panel_builds() == builds && panvar::prepared_locus_panel_precounts() == precounts,
            "the loaded-panel route executes neither the pre-count nor build_marker_panel");
    require(loaded->panel_source() == "loaded" && report.panel_source == "loaded" && report.panel_projection_bytes == 0,
            "the panel source is recorded; no panel projection is charged");
    // every field equals a freshly built panel (from a freshly reconstructed model)
    const panvar::PreparedMarkerPanel fresh = panvar::build_marker_panel(panvar::prepare_genotype_model(f.prepared.record.input), 21, 11);
    const std::string differs = panvar::prepared_marker_panel_difference(loaded->panel(), fresh);
    require(differs.empty(), "the loaded panel differs from a fresh construction at " + differs);
    require(panvar::prepared_marker_panel_difference(loaded->panel(), f.prepared.locus->panel()).empty(), "and from the prepared one");
    require(loaded->fingerprint() == f.prepared.locus->fingerprint() && loaded->digest() == f.prepared.locus->digest(), "same locus");
    // the v1 route on the v1 file: rebuilt, recorded, counted
    panvar::PreparedLocusLoadReport v1;
    const std::uint64_t builds_v1 = panvar::marker_panel_builds(), precounts_v1 = panvar::prepared_locus_panel_precounts();
    const auto rebuilt = load_bytes(f.bytes, "v1.pvpl", kDefault, &v1);
    require(rebuilt->panel_source() == "rebuilt" && v1.panel_rebuild_reason == "no panel section" &&
                panvar::marker_panel_builds() == builds_v1 + 1 && panvar::prepared_locus_panel_precounts() == precounts_v1 + 1,
            "a v1 index rebuilds through the existing path, recorded and counted");
    require(report.accounted_peak_bytes < v1.accounted_peak_bytes, "the loose reconstruction bound is gone from the panel route");
    // an unsupported extension version is not understood: rebuilt (verified by the stored digest), reason recorded,
    // and kept verbatim by a rewrite
    std::string future = p.bytes;
    const auto [at, len] = section(future, 16);
    (void)len;
    put_u32(future, at, panvar::kPanelExtensionVersion + 1);
    panvar::reseal_prepared_locus(future);
    panvar::PreparedLocusLoadReport fr;
    const auto future_loaded = load_bytes(future, "future.pvpl", kDefault, &fr);
    require(future_loaded->panel_source() == "rebuilt" && fr.panel_rebuild_reason.find("panel_extension_version 2") != std::string::npos,
            "an unsupported panel extension version is rebuilt with its reason: " + fr.panel_rebuild_reason);
    require(panvar::encode_prepared_locus(decode(future)) == future, "an unsupported section is preserved by a rewrite");
    std::string algo = p.bytes;
    put_u32(algo, at + 4, panvar::kMarkerPanelAlgorithmVersion + 1);
    panvar::reseal_prepared_locus(algo);
    panvar::PreparedLocusLoadReport ar;
    (void)load_bytes(algo, "algo.pvpl", kDefault, &ar);
    require(ar.panel_source == "rebuilt" && ar.panel_rebuild_reason.find("marker_panel_algorithm_version 2") != std::string::npos,
            "an unsupported panel-algorithm version is rebuilt with its reason");
    remove_exact(path("panel.pvpl"));
    remove_exact(path("panel_rewrite.pvpl"));
    std::cout << "  ok: panel extension (" << p.bytes.size() - f.bytes.size() << " more bytes): independent preparations and "
                 "rewrite byte-identical; loaded route ran neither pre-count nor build; every field equals a fresh build; "
                 "accounted peak " << report.accounted_peak_bytes << " vs v1 " << v1.accounted_peak_bytes
              << " B; v1 file rebuilt (recorded, counted); unsupported versions rebuilt with reasons and preserved\n";
}

std::uint64_t revcomp_code(std::uint64_t code, std::size_t k) {
    std::uint64_t out = 0;
    for (std::size_t i = 0; i < k; ++i) {
        out = (out << 2) | (3 - (code & 3));
        code >>= 2;
    }
    return out;
}

void test_panel_malformed(const PanelFixture& p) {
    const panvar::PreparedLocusRecord& base = p.prepared.record;
    const std::size_t k = 21, s = 11;
    struct Case {
        std::string what;
        std::function<void(panvar::PreparedMarkerPanel&)> change;
    };
    // locate a block with markers and an allele list with at least two entries
    std::size_t mb = 0;
    while (base.panel->blocks[mb].marker_slots.size() < 2) ++mb;
    std::size_t ma = 0;
    while (base.panel->blocks[mb].allele_markers[ma].size() < 2) ++ma;
    const auto find_code = [&](bool canonical, bool syncmer) {
        for (std::uint64_t c = base.panel->codes.back() + 1;; ++c) {
            if ((c <= revcomp_code(c, k)) == canonical && panvar::is_closed_syncmer(c, k, s) == syncmer) return c;
        }
    };
    const std::vector<Case> structural = {
        {"codes not increasing", [](auto& q) { std::swap(q.codes[0], q.codes[1]); }},
        {"code beyond 4^k", [&](auto& q) { q.codes.back() = 1ull << (2 * k); }},
        {"non-canonical code", [&](auto& q) { q.codes.back() = find_code(false, true); }},
        {"code that is not a closed syncmer", [&](auto& q) { q.codes.back() = find_code(true, false); }},
        {"panel kmer differs", [](auto& q) { q.kmer_size = 23; }},
        {"block count", [](auto& q) { q.blocks.pop_back(); }},
        {"block id", [](auto& q) { q.blocks[0].block_id += "x"; }},
        {"allele count", [](auto& q) { q.blocks[0].allele_count += 1; }},
        {"allele lists", [&](auto& q) { q.blocks[mb].allele_markers.pop_back(); }},
        {"fixed multiplicity count", [&](auto& q) { q.blocks[mb].fixed_haploid_multiplicity.pop_back(); }},
        {"marker slot out of range", [&](auto& q) { q.blocks[mb].marker_slots.back() = static_cast<std::uint32_t>(q.codes.size()); }},
        {"marker slots not increasing", [&](auto& q) { std::swap(q.blocks[mb].marker_slots[0], q.blocks[mb].marker_slots[1]); }},
        {"allele marker index out of range", [&](auto& q) {
             q.blocks[mb].allele_markers[ma].back().first = static_cast<std::uint32_t>(q.blocks[mb].marker_slots.size());
         }},
        {"allele markers not increasing", [&](auto& q) { std::swap(q.blocks[mb].allele_markers[ma][0], q.blocks[mb].allele_markers[ma][1]); }},
        {"zero copies", [&](auto& q) { q.blocks[mb].allele_markers[ma][0].second = 0; }},
        {"a marker no allele carries", [&](auto& q) {
             const std::uint32_t local = q.blocks[mb].allele_markers[ma][0].first;
             for (auto& allele : q.blocks[mb].allele_markers) {
                 allele.erase(std::remove_if(allele.begin(), allele.end(), [&](const auto& e) { return e.first == local; }), allele.end());
             }
         }},
        {"equal copies in every allele", [&](auto& q) {
             const std::uint32_t local = q.blocks[mb].allele_markers[ma][0].first;
             for (auto& allele : q.blocks[mb].allele_markers) {
                 auto it = std::lower_bound(allele.begin(), allele.end(), std::make_pair(local, 0u));
                 if (it != allele.end() && it->first == local) it->second = 1;
                 else allele.insert(it, {local, 1u});
             }
         }},
        {"anchor sizes", [](auto& q) { q.anchor_haploid_multiplicity.pop_back(); }},
        {"anchor out of range", [](auto& q) { q.anchor_slots.back() = static_cast<std::uint32_t>(q.codes.size()); }},
        {"anchors not increasing", [](auto& q) { std::swap(q.anchor_slots[0], q.anchor_slots[1]); }},
        {"anchor with zero copies", [](auto& q) { q.anchor_haploid_multiplicity[0] = 0; }},
        {"anchor that is also a block marker", [&](auto& q) { q.anchor_slots[0] = q.blocks[mb].marker_slots[0]; }},
        {"a code neither anchor nor marker", [](auto& q) { q.anchor_slots.pop_back(); q.anchor_haploid_multiplicity.pop_back(); }},
        {"total alleles", [](auto& q) { q.total_alleles += 1; }},
        {"scored alleles", [](auto& q) { q.scored_alleles -= 1; }},
    };
    const std::map<std::string, std::string> expected = {
        {"codes not increasing", "not strictly increasing"},
        {"code beyond 4^k", "exceeds 4^k"},
        {"non-canonical code", "not canonical"},
        {"code that is not a closed syncmer", "closed syncmer"},
        {"panel kmer differs", "kmer_size or syncmer_s differs"},
        {"block count", "one inventory per model block"},
        {"block id", "id differs from the model"},
        {"allele count", ": allele count"},
        {"allele lists", "allele marker lists"},
        {"fixed multiplicity count", "fixed multiplicities"},
        {"marker slot out of range", "marker slot is out of range"},
        {"marker slots not increasing", "marker slots are not increasing"},
        {"allele marker index out of range", "allele marker index is out of range"},
        {"allele markers not increasing", "allele markers are not increasing"},
        {"zero copies", "allele marker with zero copies"},
        {"a marker no allele carries", "no allele carries"},
        {"equal copies in every allele", "equal copies in every allele"},
        {"anchor sizes", "anchor slots and multiplicities"},
        {"anchor out of range", "anchor slot is out of range"},
        {"anchors not increasing", "anchor slots are not increasing"},
        {"anchor with zero copies", "anchor with zero copies"},
        {"anchor that is also a block marker", "both an anchor and a block marker"},
        {"a code neither anchor nor marker", "neither an anchor"},
        {"total alleles", "total or scored"},
        {"scored alleles", "total or scored"},
    };
    for (const Case& c : structural) {
        panvar::PreparedLocusRecord r = base;
        c.change(*r.panel);
        const std::string bytes = panvar::encode_prepared_locus(r);
        const std::string message = expect(Kind::invalid_structure, [&] { decode(bytes); }, "panel: " + c.what);
        require(message.find(expected.at(c.what)) != std::string::npos, "panel: " + c.what + " tripped another rule: " + message);
        expect(Kind::invalid_structure, [&] { load_bytes(bytes, "malformed.pvpl"); }, "panel (load, never rebuilt): " + c.what);
    }
    // fields no structural rule constrains: caught by the stored complete digest before use, never rebuilt
    const std::vector<Case> digest_only = {
        {"dropped multi-block count", [](auto& q) { q.dropped_multiblock_markers += 1; }},
        {"a fixed multiplicity value", [&](auto& q) { q.blocks[mb].fixed_haploid_multiplicity[0] += 1; }},
        {"an anchor multiplicity value", [](auto& q) { q.anchor_haploid_multiplicity[0] += 1; }},
        {"an allele's copy count", [&](auto& q) { q.blocks[mb].allele_markers[ma][0].second += 7; }},
    };
    for (const Case& c : digest_only) {
        panvar::PreparedLocusRecord r = base;
        c.change(*r.panel);
        const std::string bytes = panvar::encode_prepared_locus(r);
        (void)decode(bytes);
        expect(Kind::reconstruction_mismatch, [&] { load_bytes(bytes, "digest.pvpl"); }, "panel digest: " + c.what);
    }
    // byte-level: a different model dependency; a huge count; the section flagged required
    const auto [at, len] = section(p.bytes, 16);
    (void)len;
    {
        std::string x = p.bytes;
        x[at + 8 + 8] = x[at + 8 + 8] == 'a' ? 'b' : 'a';   // first hex digit of the model dependency
        panvar::reseal_prepared_locus(x);
        expect(Kind::invalid_structure, [&] { decode(x); }, "panel model dependency");
        std::string y = p.bytes;
        put_u64(y, at + 8 + 8 + 64 + 16, 1ull << 40);   // the codes count
        panvar::reseal_prepared_locus(y);
        expect(Kind::limit_exceeded, [&] { decode(y); }, "panel codes count beyond max_elements");
        std::string z = p.bytes;
        put_u32(z, at - 12, 1);
        panvar::reseal_prepared_locus(z);
        expect(Kind::invalid_section, [&] { decode(z); }, "section 16 flagged required");
    }
    std::cout << "  ok: " << structural.size() << " malformed panels refused (decode and load; never rebuilt); "
              << digest_only.size() << " digest-only changes refused before use; model dependency, count and required-flag damage refused\n";
}

void test_panel_accounting(const PanelFixture& p) {
    panvar::PreparedLocusLoadReport report;
    (void)load_bytes(p.bytes, "pacct.pvpl", kDefault, &report);
    std::uint64_t running = 0;
    std::size_t phases = 0;
    bool panel_phase = false;
    for (const auto& [label, peak] : report.phase_peaks) {
        if (peak <= running) continue;
        running = peak;
        panvar::PreparedLocusLoadLimits tight = kDefault;
        tight.max_accounted_bytes = peak - 1;
        const std::string message = expect(Kind::limit_exceeded, [&] { load_bytes(p.bytes, "pacct.pvpl", tight); }, "panel limit at " + label);
        require(message.find(label) != std::string::npos, "the refusal names its phase: " + label);
        panel_phase |= label.rfind("panel", 0) == 0;
        ++phases;
    }
    require(panel_phase, "the decoded panel's own charges are enforced phases");
    for (const auto& [label, peak] : report.phase_peaks) {
        (void)peak;
        require(label.find("bound") == std::string::npos && label.find("pre-count") == std::string::npos,
                "no reconstruction bound on the panel route: " + label);
    }
    panvar::PreparedLocusRecord r = p.prepared.record;
    r.panel->blocks[0].allele_markers.resize(400000);   // empty lists: 8 bytes encoded each, 24 and up decoded
    const std::string bytes = panvar::encode_prepared_locus(r);
    panvar::PreparedLocusLoadLimits tight = kDefault;
    tight.max_accounted_bytes = 2 * bytes.size();
    const std::string message = expect(Kind::limit_exceeded, [&] { decode(bytes, tight); }, "panel container amplification");
    require(message.find("panel alleles") != std::string::npos, "caught while decoding the panel: " + message);
    std::cout << "  ok: panel route accounting peak " << report.accounted_peak_bytes << " B (decoded panel "
              << report.panel_decoded_accounted_bytes << " B; no bound); refused one byte below each of " << phases
              << " phase maxima; allocation amplification inside the panel refused\n";
}

void test_hardening(const Fixture& f, const PanelFixture& p) {
    // marker settings validated before any pre-count or reconstruction: a typed refusal, v1 and panel files alike
    const std::vector<std::pair<std::size_t, std::size_t>> settings = {{32, 11}, {21, 0}, {11, 21}, {0, 0}};
    for (const auto& [k, s] : settings) {
        for (const panvar::PreparedLocusRecord* base : {&f.prepared.record, &p.prepared.record}) {
            panvar::PreparedLocusRecord r = *base;
            r.key.kmer_size = k;
            r.key.syncmer_s = s;
            if (r.panel) {
                r.panel->kmer_size = k;
                r.panel->syncmer_s = s;
            }
            const std::string bytes = panvar::encode_prepared_locus(r);
            expect(Kind::invalid_structure, [&] { decode(bytes); }, "marker settings k=" + std::to_string(k) + " s=" + std::to_string(s));
            expect(Kind::invalid_structure, [&] { load_bytes(bytes, "settings.pvpl"); }, "marker settings (load)");
        }
    }
    // post-publication failures: the target exists, is complete, is left in place, and the error says so
    const std::string target = path("post.pvpl");
    const auto temporaries = [&] {
        std::vector<std::string> out;
        for (const auto& entry : std::filesystem::directory_iterator(g_dir)) {
            if (entry.path().filename().string().rfind("post.pvpl.tmp.", 0) == 0) out.push_back(entry.path().string());
        }
        return out;
    };
    for (const int stage : {3, 4, 5}) {
        remove_exact(target);
        for (const std::string& stale : temporaries()) remove_exact(stale);
        const std::string message = expect(Kind::post_publication, [&] { panvar::publish_prepared_locus_stage(f.bytes, target, stage); },
                                           "post-publication stage " + std::to_string(stage));
        require(message.find("was published") != std::string::npos && message.find("left in place") != std::string::npos,
                "the error states the target was published: " + message);
        require(slurp(target) == f.bytes, "the published target is complete and untouched");
        require((stage == 3) == (temporaries().size() == 1), "a failed unlink leaves its temporary, and only then");
        require(panvar::load_prepared_locus(target, kDefault)->fingerprint() == f.prepared.locus->fingerprint(), "and it loads");
    }
    remove_exact(target);
    for (const std::string& stale : temporaries()) remove_exact(stale);
    std::cout << "  ok: invalid marker settings refused (typed) before reconstruction in v1 and panel files; post-publication "
                 "unlink, directory-open and directory-fsync failures reported with the valid target left in place\n";
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: genotype_prepared_locus_file_test SCRATCH_DIR\n";
        return 2;
    }
    g_dir = argv[1];
    std::filesystem::create_directories(g_dir);
    try {
        const Fixture f = build_fixture();
        std::cout << "  ok: fixture built by the real builder: curated partial template (in coverage, absent from the model), "
                     "mixed orientation, alias walk, the builder's empty-interior source\n";
        test_determinism_round_trip_rewrite(f);
        test_empty_sequence_allele(f);
        test_versions_and_provenance(f);
        test_container_corruption(f);
        test_semantic_bytes(f);
        test_semantic_records(f);
        test_accounting(f);
        test_publication(f);
        const PanelFixture p = build_panel_fixture(f);
        test_panel_round_trip(f, p);
        test_panel_malformed(p);
        test_panel_accounting(p);
        test_hardening(f, p);
        std::cout << "prepared locus file: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "prepared locus file: FAIL: " << error.what() << '\n';
        return 1;
    }
}
