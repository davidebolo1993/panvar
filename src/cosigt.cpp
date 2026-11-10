#include "panvar/cosigt.hpp"
#include "panvar/cli_utils.hpp"

#include <algorithm>
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
        if (pos == std::string::npos) { out.push_back(s.substr(start)); break; }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

} // namespace


// Which contig of a slot stands for it where `observed` says which contigs were seen: exactly one. None
// seen, or two (overlapping contigs, which would count the same sequence twice), leaves the slot -- and so
// the sample -- unknown there, even where the two would agree; call --diploid instead keeps agreeing
// overlapping contigs, since it compares their values.
const std::string* slot_contig(const std::vector<std::string>& slot,
                               const std::function<bool(const std::string&)>& observed) {
    const std::string* hit = nullptr;
    for (const std::string& c : slot) {
        if (!observed(c)) continue;
        if (hit != nullptr && *hit != c) return nullptr;
        hit = &c;
    }
    return hit;
}

// Parse a cosigt sample->haplotype table: "sample<TAB>hap1<TAB>hap2..." (header row tolerated), a
// haplotype in several contigs written "contig1,contig2" in its column. Shared by every module taking
// --samples, so a sample's haplotypes are read one way everywhere.
CosigtTable read_cosigt_table(const std::string& path,
                              const std::string& module,
                              const std::unordered_set<std::string>* known_haplotypes,
                              const char* known_label,
                              bool quiet) {
    const std::string where = module + " --samples: ";
    std::ifstream sin(path);
    if (!sin) throw std::runtime_error("Failed to open --samples file: " + path);
    // Read as plain text: a gzipped table would otherwise be read as bytes, into sample ids of noise.
    if (sin.peek() == 0x1f) {
        char magic[2] = {0, 0};
        sin.read(magic, 2);
        if (static_cast<unsigned char>(magic[1]) == 0x8b)
            throw std::runtime_error(where + path + " is gzip-compressed; give it uncompressed");
        sin.clear();
        sin.seekg(0);
    }
    CosigtTable t;
    std::unordered_set<std::string> seen_samples;
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
        t.sample_order.push_back(sample);
        std::vector<std::vector<std::string>>& slots = t.slots[sample];
        for (std::size_t fi = 1; fi < fields.size(); ++fi) {
            std::vector<std::string> slot;
            for (const std::string& hap : split_on(fields[fi], ',')) {
                if (hap.empty()) continue;
                // A name that matches nothing would contribute nothing and look like a sample with no
                // signal, so it is rejected.
                if (known_haplotypes != nullptr && !known_haplotypes->count(hap))
                    throw std::runtime_error(where + "sample '" + sample + "' names haplotype '" +
                                             hap + "', which is not " + known_label);
                if (std::find(slot.begin(), slot.end(), hap) != slot.end()) continue;
                slot.push_back(hap);
                // A haplotype in two slots of one sample is a homozygous genotype, counted twice. One in
                // several samples is allowed, as in a cohort drawn from a panel with replacement.
                t.path_to_samples[hap].push_back(sample);
            }
            if (!slot.empty()) slots.push_back(std::move(slot));
        }
        if (slots.empty())
            throw std::runtime_error(where + "sample '" + sample + "' has no haplotype "
                                     "assignment; a row with no haplotypes cannot be genotyped");
    }
    if (t.sample_order.empty())
        throw std::runtime_error(where + path + " lists no sample: a header alone gives a sample level with no "
                                 "columns");
    // Ploidy is not assumed: a sample's dosage is the sum over its own haplotypes, so any number works.
    // Mixed ploidy is noted rather than refused.
    std::size_t lo = 0, hi = 0;
    for (const auto& [sm, sl] : t.slots) {
        if (lo == 0 || sl.size() < lo) lo = sl.size();
        if (sl.size() > hi) hi = sl.size();
    }
    if (lo != hi && !quiet)
        std::cerr << cli::log_prefix(module) << "note: --samples has mixed ploidy (" << lo << ".." << hi
                  << " haplotypes per sample); dosage is the sum over each sample's own haplotypes\n";
    return t;
}

} // namespace panvar
