// Private k-mers of locus indexes found in a FASTA outside listed contigs: what a background holding
// only some chromosomes misses. Per locus: private k-mers and how many occur elsewhere; with
// PSCAN_DUMP=<prefix>, their codes to <prefix>.<locus>.found (for dev/relabel_index.cpp).
//   c++ -std=c++20 -O2 -Iinclude dev/private_elsewhere.cpp -lz -o private_elsewhere
//   private_elsewhere <fasta.gz> <skip contigs,comma> <locus>=<private k-mers, one per line> ...
#include "panvar/kmer.hpp"
#include <zlib.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
using namespace panvar;
int main(int argc, char** argv) {
    // pscan <fasta.gz> <skip contigs,comma> <locus=P.txt>...
    std::set<std::string> skip;
    { std::string s = argv[2]; size_t p = 0; while (p <= s.size()) { size_t q = s.find(',', p); if (q == std::string::npos) q = s.size(); skip.insert(s.substr(p, q - p)); p = q + 1; } }
    std::vector<std::string> names;
    std::unordered_map<std::uint64_t, std::uint32_t> loc;
    std::vector<std::size_t> total;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i]; auto eq = a.find('=');
        names.push_back(a.substr(0, eq)); total.push_back(0);
        std::ifstream f(a.substr(eq + 1)); std::string km;
        while (f >> km) { kmer::for_each_canonical_kmer(km.data(), km.size(), 31, [&](std::uint64_t c, std::size_t, bool) { loc[c] |= 1u << (i - 3); }); ++total.back(); }
    }
    std::unordered_map<std::uint64_t, std::uint32_t> seen;
    gzFile g = gzopen(argv[1], "rb"); gzbuffer(g, 1 << 20);
    std::vector<char> buf(1 << 16); std::string seq, name; std::size_t scanned = 0;
    auto flush = [&] {
        if (!name.empty() && !skip.count(name)) {
            kmer::for_each_canonical_kmer(seq.data(), seq.size(), 31, [&](std::uint64_t c, std::size_t, bool) { auto it = loc.find(c); if (it != loc.end()) ++seen[c]; });
            scanned += seq.size();
        }
        seq.clear();
    };
    while (gzgets(g, buf.data(), buf.size())) {
        if (buf[0] == '>') { flush(); name = std::string(buf.data() + 1); name = name.substr(0, name.find_first_of(" \t\n")); continue; }
        std::size_t n = std::strlen(buf.data()); while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) --n; seq.append(buf.data(), n);
    }
    flush();
    std::printf("scanned %.2f Gb outside the background\n", scanned / 1e9);
    if (const char* dump = std::getenv("PSCAN_DUMP")) {
        for (std::size_t i = 0; i < names.size(); ++i) {
            std::ofstream o(std::string(dump) + "." + names[i] + ".found", std::ios::app);
            for (auto& [c, n] : seen) if (loc[c] >> i & 1) o << c << "\n";
        }
    }
    std::printf("locus\tprivate\tfound_elsewhere\tshare\n");
    for (std::size_t i = 0; i < names.size(); ++i) {
        std::size_t f = 0; for (auto& [c, n] : seen) if (loc[c] >> i & 1) ++f;
        std::printf("%s\t%zu\t%zu\t%.4f\n", names[i].c_str(), total[i], f, total[i] ? double(f) / total[i] : 0.0);
    }
}
