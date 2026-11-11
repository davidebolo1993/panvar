#include "panvar/recruit_reads.hpp"

#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace panvar::recruit {

namespace {

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint8_t kLabelP = 0, kLabelS = 1, kLabelR = 2, kLabelD = 3, kLabelN = 4, kMultiBit = 8;

struct Fragment {
    std::string name, seq1, seq2;
};

// One locus, as reads see it: each k-mer's label (P/S/R, D for decoys), with a bit for multi-copy P.
struct LocusTable {
    std::string name;
    std::unordered_map<std::uint64_t, std::uint8_t> label;
};

struct Row {
    std::string fragment;
    std::uint32_t locus = 0;
    Evidence ev;
    const char* fate = "";
    const char* reason = "";
};

// Labels of one read's k-mers against a locus. Private and decoy k-mers count as independent sites only
// k bases apart (greedily, in read order), so one differing base, whose k-mers may be split by a shared
// one, is one site; codes counted in the other mate (`seen`) are not counted again. Runs of unlabelled
// k-mers are counted; a base other than ACGT ends runs and is not itself evidence.
void label_read(const std::string& s, const LocusTable& t, Evidence& ev, std::vector<std::uint64_t>& seen, bool second) {
    ev.bases += static_cast<std::uint32_t>(s.size());
    if (s.size() >= kK) ev.windows += static_cast<std::uint32_t>(s.size() - kK + 1);
    std::size_t last_p = 0, last_d = 0, last_start = 0;
    bool have_p = false, have_d = false, in_n = false;
    auto counted_before = [&](std::uint64_t c) { return second && std::binary_search(seen.begin(), seen.end(), c); };
    std::vector<std::uint64_t> mine;
    kmer::for_each_canonical_kmer(s.data(), s.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
        ++ev.kmers;
        const auto it = t.label.find(c);
        int lab = kLabelN;
        if (it != t.label.end()) {
            lab = it->second & 7;
            if (lab == kLabelP && (it->second & kMultiBit)) ++ev.multi;
        }
        if (in_n && (lab != kLabelN || start != last_start + 1)) in_n = false;
        switch (lab) {
            case kLabelP:
                ++ev.n_p;
                if (!counted_before(c) && (!have_p || start >= last_p + kK)) { ++ev.sites_p; last_p = start; have_p = true; }
                mine.push_back(c);
                break;
            case kLabelD:
                ++ev.n_d;
                if (!counted_before(c) && (!have_d || start >= last_d + kK)) { ++ev.sites_d; last_d = start; have_d = true; }
                mine.push_back(c);
                break;
            case kLabelS: ++ev.n_s; break;
            case kLabelR: ++ev.n_r; break;
            default:
                ++ev.n_n;
                if (!in_n) { ++ev.runs_n; in_n = true; }
        }
        last_start = start;
    });
    if (!second) {
        seen.swap(mine);
        std::sort(seen.begin(), seen.end());
    }
}

// The most mismatch runs a locus read of this many callable k-mers should show: an upper Poisson bound on
// sequencing errors plus about one new allele per kilobase (normal approximation, 3.1 standard deviations).
double tolerated_runs(std::uint32_t callable, double e) {
    const double lambda = (e + 0.001) * static_cast<double>(callable);
    return lambda + 3.1 * std::sqrt(lambda) + 2.0;
}

void decide(Row& r, const ReadOptions& o) {
    const Evidence& v = r.ev;
    const bool uncallable = v.windows > 0 && v.kmers * 2 < v.windows;
    const bool divergent = static_cast<double>(v.runs_n) > tolerated_runs(v.kmers, o.error_rate);
    bool certain = false;
    if (!uncallable && !divergent) {
        switch (o.method) {
            case Method::Runs:
                certain = v.sites_p >= o.min_sites && v.sites_d <= o.decoy_share * v.sites_p &&
                          v.sites_p >= o.min_density * static_cast<double>(v.kmers) / static_cast<double>(kK);
                break;
            case Method::Kmers: certain = v.n_d == 0 && v.n_p >= std::max<std::uint32_t>(3, o.min_sites); break;
            case Method::Fraction:
                certain = v.kmers > 0 && v.n_p >= 1 && static_cast<double>(v.n_p + v.n_s + v.n_r) / v.kmers >= o.min_fraction;
                break;
        }
    }
    if (certain) { r.fate = "certain"; r.reason = "."; return; }
    if (v.sites_p == 0 && v.sites_d >= 2) { r.fate = "unrecruited"; r.reason = "homolog"; return; }
    if (v.sites_p >= 1 && v.sites_d >= 1) { r.fate = "ambiguous"; r.reason = "conflict"; return; }
    if (uncallable && v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "uncallable"; return; }
    if (divergent && v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "divergent"; return; }
    if (v.sites_p >= o.min_sites) { r.fate = "ambiguous"; r.reason = "sparse"; return; }
    if (v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "single_site"; return; }
    if (v.n_s >= 1) { r.fate = "ambiguous"; r.reason = "shared_only"; return; }
    r.fate = "unrecruited";
    r.reason = "no_evidence";
}

bool read_fastq(GzLineReader& in, std::string& name, std::string& seq) {
    std::string h, plus, qual;
    while (in.getline(h)) {
        if (h.empty()) continue;
        if (h[0] != '@') throw std::runtime_error("recruit: not FASTQ (a record starts with '" + h.substr(0, 1) + "')");
        if (!in.getline(seq) || !in.getline(plus) || !in.getline(qual)) throw std::runtime_error("recruit: a truncated FASTQ record");
        if (!seq.empty() && seq.back() == '\r') seq.pop_back();
        if (!qual.empty() && qual.back() == '\r') qual.pop_back();
        if (plus.empty() || plus[0] != '+' || qual.size() != seq.size())
            throw std::runtime_error("recruit: a malformed FASTQ record at " + h.substr(0, 60) + " (multi-line FASTQ is not read)");
        const auto sp = h.find_first_of(" \t");
        name = h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
        if (name.size() > 2 && name[name.size() - 2] == '/' && (name.back() == '1' || name.back() == '2')) name.resize(name.size() - 2);
        return true;
    }
    return false;
}

}  // namespace

void Evidence::add(const Evidence& o) {
    sites_p += o.sites_p; sites_d += o.sites_d; runs_n += o.runs_n;
    n_p += o.n_p; n_s += o.n_s; n_r += o.n_r; n_d += o.n_d; n_n += o.n_n;
    kmers += o.kmers; windows += o.windows; bases += o.bases; multi += o.multi;
}

ReadStats recruit_reads(const ReadOptions& o) {
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    if (o.indexes.empty()) throw std::runtime_error("recruit: no locus index");
    if (o.indexes.size() > 64) throw std::runtime_error("recruit: at most 64 loci in one run");
    std::vector<LocusTable> loci(o.indexes.size());
    std::unordered_map<std::uint64_t, std::uint64_t> screen;   // syncmer -> loci it belongs to
    for (std::size_t li = 0; li < o.indexes.size(); ++li) {
        const LocusIndex x = load(o.indexes[li]);
        LocusTable& t = loci[li];
        t.name = x.name;
        t.label.reserve(x.codes.size() + x.decoy.size());
        for (std::size_t i = 0; i < x.codes.size(); ++i)
            t.label[x.codes[i]] = static_cast<std::uint8_t>(x.cls[i] | (x.cls[i] == kPrivate && x.mult[i] > 1 ? kMultiBit : 0));
        for (const std::uint64_t d : x.decoy) t.label.emplace(d, kLabelD);
        for (const std::uint64_t s : x.screen) screen[s] |= 1ULL << li;
        log(x.name + ": " + std::to_string(x.codes.size()) + " k-mers, " + std::to_string(x.decoy.size()) + " decoys, " +
            std::to_string(x.screen.size()) + " screen syncmers");
    }

    GzLineReader in1(o.reads1);
    if (!in1.ok()) throw std::runtime_error("recruit: cannot read " + o.reads1);
    std::unique_ptr<GzLineReader> in2;
    if (!o.reads2.empty()) {
        in2 = std::make_unique<GzLineReader>(o.reads2);
        if (!in2->ok()) throw std::runtime_error("recruit: cannot read " + o.reads2);
    }
    std::FILE* out = std::fopen(o.out.c_str(), "w");
    if (out == nullptr) throw std::runtime_error("recruit: cannot write " + o.out);
    std::fprintf(out, "fragment\tlocus\tfate\treason\tsites_p\tsites_d\truns_n\tn_p\tn_s\tn_r\tn_d\tn_n\tkmers\twindows\tbases\tmulti_p\n");

    ReadStats st;
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    constexpr std::size_t kBatch = 20000;
    std::vector<Fragment> batch;
    std::vector<std::vector<Row>> rows;
    while (true) {
        batch.clear();
        Fragment f;
        std::string n2;
        while (batch.size() < kBatch && read_fastq(in1, f.name, f.seq1)) {
            if (in2) {
                if (!read_fastq(*in2, n2, f.seq2)) throw std::runtime_error("recruit: " + o.reads2 + " has fewer reads than " + o.reads1);
                if (n2 != f.name) throw std::runtime_error("recruit: mates out of step: " + f.name + " and " + n2);
            }
            batch.push_back(f);
        }
        if (batch.empty()) break;
        rows.assign(batch.size(), {});
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                for (std::size_t i; (i = next++) < batch.size();) {
                    const Fragment& fr = batch[i];
                    std::uint64_t hit = 0;
                    for (const std::string* s : {&fr.seq1, &fr.seq2})
                        kmer::for_each_closed_syncmer(s->data(), s->size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
                            const auto it = screen.find(c);
                            if (it != screen.end()) hit |= it->second;
                        });
                    for (std::uint32_t li = 0; hit != 0; ++li, hit >>= 1) {
                        if (!(hit & 1)) continue;
                        Row r;
                        r.fragment = fr.name;
                        r.locus = li;
                        std::vector<std::uint64_t> seen;
                        label_read(fr.seq1, loci[li], r.ev, seen, false);
                        if (!fr.seq2.empty()) label_read(fr.seq2, loci[li], r.ev, seen, true);
                        decide(r, o);
                        rows[i].push_back(r);
                    }
                }
            });
        for (auto& t : pool) t.join();
        for (std::size_t i = 0; i < batch.size(); ++i) {
            ++st.fragments;
            if (rows[i].empty()) continue;
            ++st.screened;
            for (const Row& r : rows[i]) {
                const Evidence& v = r.ev;
                std::fprintf(out, "%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", r.fragment.c_str(), loci[r.locus].name.c_str(), r.fate,
                             r.reason, v.sites_p, v.sites_d, v.runs_n, v.n_p, v.n_s, v.n_r, v.n_d, v.n_n, v.kmers, v.windows, v.bases, v.multi);
                if (r.fate[0] == 'c') ++st.certain;
                else if (r.fate[0] == 'a') ++st.ambiguous;
                else ++st.unrecruited;
            }
        }
    }
    if (in2) {
        std::string n2, s2;
        if (read_fastq(*in2, n2, s2)) throw std::runtime_error("recruit: " + o.reads2 + " has more reads than " + o.reads1);
    }
    if (std::fclose(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);
    return st;
}

}  // namespace panvar::recruit
