// Relabels the private k-mers listed in a file (codes, one per line) as shared, as a background that also
// held those sequences would, and saves the index under a new name: an emulated whole-genome index.
//   c++ -std=c++20 -O2 -Iinclude dev/relabel_index.cpp build/libpanvarlib.a build/libsshash.a build/libedlib.a -lz -o relabel_index
//   relabel_index <in.plx> <codes> <out.plx>
#include "panvar/locus_index.hpp"
#include <cstdio>
#include <fstream>
#include <unordered_set>
using namespace panvar::recruit;
int main(int argc, char** argv) {
    LocusIndex x = load(argv[1]);
    std::unordered_set<std::uint64_t> found;
    std::ifstream f(argv[2]);
    for (std::uint64_t c; f >> c;) found.insert(c);
    std::size_t moved = 0;
    for (std::size_t i = 0; i < x.codes.size(); ++i)
        if (x.cls[i] == kPrivate && found.count(x.codes[i])) { x.cls[i] = kShared; if (x.outside[i] == 0) x.outside[i] = 1; ++moved; }
    save(x, argv[3]);
    std::printf("%s: %zu private k-mers relabelled shared\n", argv[1], moved);
}
