// Edit distances of sequence pairs, for the transcript checks: reads "id<TAB>a<TAB>b" lines and prints
// "id<TAB>distance" (global, unit costs; edlib).
//   c++ -std=c++17 -O2 -Iexternal/edlib/edlib/include dev/edit_distance.cpp build/libedlib.a -o build/edit_distance
#include <iostream>
#include <string>

#include "edlib.h"

int main() {
    std::ios::sync_with_stdio(false);
    std::string line;
    while (std::getline(std::cin, line)) {
        const std::size_t t1 = line.find('\t'), t2 = t1 == std::string::npos ? t1 : line.find('\t', t1 + 1);
        if (t2 == std::string::npos) continue;
        const std::string a = line.substr(t1 + 1, t2 - t1 - 1), b = line.substr(t2 + 1);
        EdlibAlignResult r = edlibAlign(a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()),
                                        edlibDefaultAlignConfig());
        std::cout << line.substr(0, t1) << '\t' << r.editDistance << '\n';
        edlibFreeAlignResult(r);
    }
}
