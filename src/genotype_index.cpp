#include "panvar/genotype_index.hpp"

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace panvar {
namespace {

// Bumped to 5 when marker_clumps, node_first_pos, all_kmers and fragment_len joined the serialized
// panel. An index written by an older build has no such fields, and reading one as if it did would
// desynchronise every vector after it. marker_clumps in particular is not cosmetic: absent, the
// emission falls back to span/fragment_len for its effective sample size, so an indexed run silently
// used a different -- and worse -- ESS discount than the direct one.
constexpr char kMagic[8] = {'p', 'v', 'g', 't', 'i', 'd', 'x', '5'};

template <typename T>
void put(std::ostream& o, const T& v) {
    o.write(reinterpret_cast<const char*>(&v), sizeof(T));
}
template <typename T>
void get(std::istream& i, T& v) {
    i.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!i) throw std::runtime_error("genotype index: truncated file");
}
void put_str(std::ostream& o, const std::string& s) {
    put<std::uint64_t>(o, s.size());
    o.write(s.data(), static_cast<std::streamsize>(s.size()));
}
std::string get_str(std::istream& i) {
    std::uint64_t n = 0;
    get(i, n);
    std::string s(n, '\0');
    if (n) i.read(&s[0], static_cast<std::streamsize>(n));
    if (!i) throw std::runtime_error("genotype index: truncated file");
    return s;
}
template <typename T>
void put_vec(std::ostream& o, const std::vector<T>& v) {
    put<std::uint64_t>(o, v.size());
    if (!v.empty()) o.write(reinterpret_cast<const char*>(v.data()),
                            static_cast<std::streamsize>(v.size() * sizeof(T)));
}
template <typename T>
void get_vec(std::istream& i, std::vector<T>& v) {
    std::uint64_t n = 0;
    get(i, n);
    v.resize(n);
    if (n) i.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
    if (!i) throw std::runtime_error("genotype index: truncated file");
}

} // namespace

void write_genotype_index(const std::string& path, const GenotypeIndex& index) {
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("genotype index: cannot write " + path);
    o.write(kMagic, sizeof(kMagic));
    put<std::uint64_t>(o, index.kmer_size);
    put<std::uint64_t>(o, index.syncmer_s);

    put<std::uint64_t>(o, index.haplotype_names.size());
    for (const std::string& n : index.haplotype_names) put_str(o, n);

    put<std::uint64_t>(o, index.chain.size());
    for (const Block& b : index.chain) {
        put<std::uint64_t>(o, b.index);
        put<std::uint32_t>(o, static_cast<std::uint32_t>(b.kind));
        put<std::uint64_t>(o, b.bubble_id);
        put_str(o, b.source);
        put_str(o, b.sink);
    }

    put<std::uint64_t>(o, index.blocks.size());
    for (const BlockAlleles& b : index.blocks) {
        put<std::uint64_t>(o, b.block_index);
        put<std::uint64_t>(o, b.n_alleles);
        put<std::uint64_t>(o, b.n_walk_alleles);
        put<std::int64_t>(o, b.bypass_allele);
        put<std::uint64_t>(o, b.n_traversing);
        put_vec(o, b.allele_haplotypes);
        put_vec(o, b.allele_bp);
        put<std::uint64_t>(o, b.allele_of.size());
        for (const auto& [name, ai] : b.allele_of) { put_str(o, name); put<std::uint64_t>(o, ai); }
    }

    put_vec(o, index.panel.node_codes);
    put_vec(o, index.panel.edge_keys);
    put_vec(o, index.panel.all_edge_keys);
    put<std::uint64_t>(o, index.panel.by_block.size());
    for (const auto& per_allele : index.panel.by_block) {
        put<std::uint64_t>(o, per_allele.size());
        for (const AlleleMarkerSet& m : per_allele) {
            put_vec(o, m.nodes);
            put_vec(o, m.edges);
        }
    }
    put<std::uint64_t>(o, index.panel.anchor_slots.size());
    for (const auto& v : index.panel.anchor_slots) put_vec(o, v);
    put_vec(o, index.panel.marker_clumps);
    put_vec(o, index.panel.node_first_pos);
    put<std::uint8_t>(o, index.panel.all_kmers ? 1 : 0);
    put<double>(o, index.panel.fragment_len);
    if (!o) throw std::runtime_error("genotype index: write failed for " + path);
}

GenotypeIndex read_genotype_index(const std::string& path) {
    std::ifstream i(path, std::ios::binary);
    if (!i) throw std::runtime_error("genotype index: cannot read " + path);
    char magic[sizeof(kMagic)];
    i.read(magic, sizeof(magic));
    if (!i || std::string(magic, sizeof(magic)) != std::string(kMagic, sizeof(kMagic))) {
        throw std::runtime_error("genotype index: bad magic in " + path + " (rebuild it)");
    }

    GenotypeIndex x;
    std::uint64_t n = 0;
    get(i, n); x.kmer_size = n;
    get(i, n); x.syncmer_s = n;

    get(i, n);
    x.haplotype_names.reserve(n);
    for (std::uint64_t k = 0; k < n; ++k) x.haplotype_names.push_back(get_str(i));

    get(i, n);
    x.chain.resize(n);
    for (Block& b : x.chain) {
        std::uint64_t v = 0;
        std::uint32_t kind = 0;
        get(i, v); b.index = v;
        get(i, kind); b.kind = static_cast<BlockKind>(kind);
        get(i, v); b.bubble_id = v;
        b.source = get_str(i);
        b.sink = get_str(i);
    }

    get(i, n);
    x.blocks.resize(n);
    for (BlockAlleles& b : x.blocks) {
        std::uint64_t v = 0;
        get(i, v); b.block_index = v;
        get(i, v); b.n_alleles = v;
        get(i, v); b.n_walk_alleles = v;
        std::int64_t bp = -1;
        get(i, bp); b.bypass_allele = static_cast<int>(bp);
        get(i, v); b.n_traversing = v;
        get_vec(i, b.allele_haplotypes);
        get_vec(i, b.allele_bp);
        get(i, v);
        for (std::uint64_t k = 0; k < v; ++k) {
            const std::string name = get_str(i);
            std::uint64_t ai = 0;
            get(i, ai);
            b.allele_of.emplace(name, ai);
        }
    }

    get_vec(i, x.panel.node_codes);
    get_vec(i, x.panel.edge_keys);
    get_vec(i, x.panel.all_edge_keys);
    get(i, n);
    x.panel.by_block.resize(n);
    for (auto& per_allele : x.panel.by_block) {
        std::uint64_t m = 0;
        get(i, m);
        per_allele.resize(m);
        for (AlleleMarkerSet& s : per_allele) {
            get_vec(i, s.nodes);
            get_vec(i, s.edges);
        }
    }
    get(i, n);
    x.panel.anchor_slots.resize(n);
    for (auto& v : x.panel.anchor_slots) get_vec(i, v);
    get_vec(i, x.panel.marker_clumps);
    get_vec(i, x.panel.node_first_pos);
    std::uint8_t ak = 0;
    get(i, ak); x.panel.all_kmers = ak != 0;
    get(i, x.panel.fragment_len);

    x.panel.kmer_size = x.kmer_size;
    x.panel.syncmer_s = x.syncmer_s;
    return x;
}

} // namespace panvar
