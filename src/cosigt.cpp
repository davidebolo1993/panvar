#include "panvar/cosigt.hpp"

#include <cctype>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace panvar {
namespace {

std::vector<std::string> split_on(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t pos = s.find(delim, start);
        if (pos == std::string::npos) {
            out.push_back(s.substr(start));
            break;
        }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

} // namespace

std::pair<std::unordered_map<std::string, std::vector<std::string>>, std::vector<std::string>>
read_cosigt_table(const std::string& path,
                  const std::string& module,
                  const std::unordered_set<std::string>* known_haplotypes,
                  const char* known_label) {
    const std::string where = module + " --samples: ";
    std::ifstream sin(path);
    if (!sin) throw std::runtime_error("Failed to open --samples file: " + path);
    std::unordered_map<std::string, std::vector<std::string>> path_to_samples;
    std::vector<std::string> sample_order;
    std::unordered_set<std::string> seen_samples;
    std::unordered_map<std::string, std::size_t> ploidy;
    std::string line;
    bool first = true;
    while (std::getline(sin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> fields = split_on(line, '\t');
        if (fields.empty()) continue;
        if (first) {
            first = false;
            std::string h = fields[0];
            for (char& c : h) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (h == "sample" || h == "sample.id" || h == "sample_id" || (!h.empty() && h[0] == '#')) continue;
        }
        const std::string& sample = fields[0];
        if (sample.empty())
            throw std::runtime_error(where + path + " has a row with an empty sample id");
        if (!seen_samples.insert(sample).second)
            throw std::runtime_error(where + "duplicate sample id '" + sample + "' in " + path +
                                     "; which row's haplotypes a sample gets would depend on file order");
        sample_order.push_back(sample);
        std::size_t assigned = 0;
        for (std::size_t fi = 1; fi < fields.size(); ++fi)
            for (const std::string& hap : split_on(fields[fi], ',')) {
                if (hap.empty()) continue;
                // A name that matches nothing contributes silently nothing, which looks exactly like a
                // sample with no signal. Reject it instead.
                if (known_haplotypes != nullptr && !known_haplotypes->count(hap))
                    throw std::runtime_error(where + "sample '" + sample + "' names haplotype '" +
                                             hap + "', which is not " + known_label);
                path_to_samples[hap].push_back(sample);
                ++assigned;
            }
        if (assigned == 0)
            throw std::runtime_error(where + "sample '" + sample + "' has no haplotype "
                                     "assignment; a row with no haplotypes cannot be genotyped");
        ploidy[sample] = assigned;
    }
    // Ploidy is not assumed. The dosage is a SUM over a sample's assigned haplotypes and every one of
    // them must be observed, so any fixed arity works -- but a table that is not uniformly diploid is
    // worth saying out loud rather than silently averaging over.
    std::size_t lo = 0, hi = 0;
    for (const auto& [sm, n] : ploidy) {
        if (lo == 0 || n < lo) lo = n;
        if (n > hi) hi = n;
    }
    if (lo != hi)
        std::cerr << "[" << module << "] note: --samples has mixed ploidy (" << lo << ".." << hi
                  << " haplotypes per sample); dosage is the sum over each sample's own assignments\n";
    return {std::move(path_to_samples), std::move(sample_order)};
}

} // namespace panvar
