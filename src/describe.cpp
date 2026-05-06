#include "panvar/describe.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <zlib.h>

namespace panvar {
namespace {

struct GeneInterval {
    std::size_t start0 = 0;
    std::size_t end0 = 0;
    std::string gene_name;
};

struct GeneModel {
    std::string chrom;
    std::string gene_name;
    std::size_t start0 = std::numeric_limits<std::size_t>::max();
    std::size_t end0 = 0;
    std::vector<std::pair<std::size_t, std::size_t>> exons;
};

struct GeneFeatureIndex {
    std::unordered_map<std::string, std::vector<GeneInterval>> gene;
    std::unordered_map<std::string, std::vector<GeneInterval>> exon;
    std::unordered_map<std::string, std::vector<GeneInterval>> cds;
    std::unordered_map<std::string, std::vector<GeneInterval>> utr;
    std::unordered_map<std::string, std::vector<GeneInterval>> intron;

    bool empty() const {
        return gene.empty() && exon.empty() && cds.empty() && utr.empty() && intron.empty();
    }
};

struct EventGeneHits {
    bool gene = false;
    bool exon = false;
    bool cds = false;
    bool utr = false;
    bool intron = false;
    std::vector<std::string> genes;
};

struct EventRecord {
    std::size_t bubble_id = 0;
    std::string record_id;
    std::string chrom;
    std::size_t pos1 = 0;
    std::size_t end1 = 0;
    std::string svtype;
    std::string event_type;
    std::string ins_subtype;
    long long svlen = 0;
    bool has_dup = false;
    std::vector<std::size_t> carrier_indices;
    EventGeneHits gene_hits;
};

struct BubbleData {
    std::vector<EventRecord> events;
};

struct HapFeatures {
    std::size_t n_events = 0;

    bool has_ins = false;
    bool has_del = false;
    bool has_inv = false;
    bool has_dup = false;

    std::size_t n_ins = 0;
    std::size_t n_del = 0;
    std::size_t n_inv = 0;
    std::size_t n_dup = 0;

    std::size_t n_events_gene = 0;
    std::size_t n_events_exon = 0;
    std::size_t n_events_cds = 0;
    std::size_t n_events_utr = 0;
    std::size_t n_events_intron = 0;
    std::unordered_set<std::string> genes_any;

    std::vector<std::size_t> ins_bins;
    std::vector<std::size_t> del_bins;
    std::vector<std::size_t> inv_bins;
    std::vector<std::size_t> dup_bins;
};

struct ParsedVcf {
    std::vector<std::string> samples;
    std::vector<std::size_t> bubble_order;
    std::unordered_map<std::size_t, BubbleData> bubbles;
};

std::string trim_ascii_whitespace(const std::string& text) {
    std::size_t lo = 0;
    while (lo < text.size() && std::isspace(static_cast<unsigned char>(text[lo]))) {
        ++lo;
    }
    std::size_t hi = text.size();
    while (hi > lo && std::isspace(static_cast<unsigned char>(text[hi - 1]))) {
        --hi;
    }
    return text.substr(lo, hi - lo);
}

std::vector<std::string> split_string(const std::string& text, char delim) {
    std::vector<std::string> out;
    std::string token;
    std::istringstream iss(text);
    while (std::getline(iss, token, delim)) {
        out.push_back(token);
    }
    return out;
}

std::string join_csv(const std::vector<std::string>& items) {
    std::ostringstream out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << items[i];
    }
    return out.str();
}

std::size_t parse_size_or_throw(const std::string& value, const std::string& name) {
    try {
        return static_cast<std::size_t>(std::stoull(value));
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid " + name + " value: " + value);
    }
}

long long parse_signed_or_default(const std::string& value, long long fallback) {
    if (value.empty()) {
        return fallback;
    }
    try {
        return std::stoll(value);
    } catch (const std::exception&) {
        return fallback;
    }
}

std::unordered_map<std::string, std::string> parse_info_field(const std::string& info) {
    std::unordered_map<std::string, std::string> out;
    const auto items = split_string(info, ';');
    for (const auto& raw_item : items) {
        if (raw_item.empty()) {
            continue;
        }
        const auto eq = raw_item.find('=');
        if (eq == std::string::npos) {
            out.emplace(raw_item, "");
            continue;
        }
        out.emplace(raw_item.substr(0, eq), raw_item.substr(eq + 1));
    }
    return out;
}

bool gt_has_alt_allele(const std::string& gt) {
    if (gt.empty() || gt == "." || gt == "./." || gt == ".|.") {
        return false;
    }

    std::string token;
    for (char c : gt) {
        if (c == '/' || c == '|') {
            if (!token.empty() && token != ".") {
                try {
                    if (std::stoll(token) > 0) {
                        return true;
                    }
                } catch (const std::exception&) {
                }
            }
            token.clear();
        } else {
            token.push_back(c);
        }
    }
    if (!token.empty() && token != ".") {
        try {
            return std::stoll(token) > 0;
        } catch (const std::exception&) {
            return false;
        }
    }
    return false;
}

std::string normalize_chrom(const std::string& chrom) {
    std::string out = trim_ascii_whitespace(chrom);
    if (out.size() >= 3 &&
        (out[0] == 'c' || out[0] == 'C') &&
        (out[1] == 'h' || out[1] == 'H') &&
        (out[2] == 'r' || out[2] == 'R')) {
        out[0] = 'c';
        out[1] = 'h';
        out[2] = 'r';
    }
    return out;
}

std::string gtf_attribute_value(const std::string& attributes, const std::string& key) {
    const std::string needle = key + " ";
    std::size_t pos = attributes.find(needle);
    if (pos == std::string::npos) {
        return {};
    }
    pos += needle.size();
    while (pos < attributes.size() && std::isspace(static_cast<unsigned char>(attributes[pos]))) {
        ++pos;
    }
    if (pos >= attributes.size()) {
        return {};
    }
    if (attributes[pos] == '"') {
        const std::size_t end = attributes.find('"', pos + 1);
        if (end == std::string::npos) {
            return {};
        }
        return attributes.substr(pos + 1, end - pos - 1);
    }
    const std::size_t end = attributes.find(';', pos);
    return trim_ascii_whitespace(attributes.substr(pos, end == std::string::npos ? std::string::npos : (end - pos)));
}

void sort_and_unique_intervals(std::unordered_map<std::string, std::vector<GeneInterval>>& by_chrom) {
    for (auto& [chrom, intervals] : by_chrom) {
        std::sort(intervals.begin(), intervals.end(), [](const GeneInterval& a, const GeneInterval& b) {
            if (a.start0 != b.start0) {
                return a.start0 < b.start0;
            }
            if (a.end0 != b.end0) {
                return a.end0 < b.end0;
            }
            return a.gene_name < b.gene_name;
        });
        intervals.erase(
            std::unique(intervals.begin(), intervals.end(), [](const GeneInterval& a, const GeneInterval& b) {
                return a.start0 == b.start0 && a.end0 == b.end0 && a.gene_name == b.gene_name;
            }),
            intervals.end());
    }
}

void add_interval(
    std::unordered_map<std::string, std::vector<GeneInterval>>& by_chrom,
    const std::string& chrom,
    std::size_t start0,
    std::size_t end0,
    const std::string& gene_name) {

    if (chrom.empty() || end0 <= start0) {
        return;
    }
    GeneInterval interval;
    interval.start0 = start0;
    interval.end0 = end0;
    interval.gene_name = gene_name.empty() ? std::string(".") : gene_name;
    by_chrom[chrom].push_back(std::move(interval));
}

std::vector<std::pair<std::size_t, std::size_t>> merge_intervals(
    std::vector<std::pair<std::size_t, std::size_t>> intervals) {

    if (intervals.empty()) {
        return {};
    }
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<std::size_t, std::size_t>> merged;
    merged.reserve(intervals.size());
    merged.push_back(intervals.front());
    for (std::size_t i = 1; i < intervals.size(); ++i) {
        auto& back = merged.back();
        const auto& curr = intervals[i];
        if (curr.first <= back.second) {
            back.second = std::max(back.second, curr.second);
        } else {
            merged.push_back(curr);
        }
    }
    return merged;
}

bool has_gzip_suffix(const std::string& path) {
    if (path.size() < 3) {
        return false;
    }
    const std::string suffix = path.substr(path.size() - 3);
    return suffix == ".gz" || suffix == ".GZ";
}

void iterate_text_lines(
    const std::string& path,
    const std::function<void(const std::string&)>& on_line,
    const std::string& readable_name) {

    if (has_gzip_suffix(path)) {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (gz == nullptr) {
            throw std::runtime_error("Failed to read gzipped " + readable_name + ": " + path);
        }
        std::string chunk(1 << 15, '\0');
        std::string line;
        while (gzgets(gz, chunk.data(), static_cast<int>(chunk.size())) != nullptr) {
            line.assign(chunk.c_str());
            while (!line.empty() && line.back() != '\n' && !gzeof(gz)) {
                if (gzgets(gz, chunk.data(), static_cast<int>(chunk.size())) == nullptr) {
                    break;
                }
                line.append(chunk.c_str());
            }
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            on_line(line);
        }

        int gz_err = Z_OK;
        const char* gz_msg = gzerror(gz, &gz_err);
        gzclose(gz);
        if (gz_err != Z_OK && gz_err != Z_STREAM_END) {
            const std::string err_text = (gz_msg == nullptr) ? "unknown zlib error" : std::string(gz_msg);
            throw std::runtime_error("Failed while reading gzipped " + readable_name + " '" + path + "': " + err_text);
        }
        return;
    }

    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to read " + readable_name + ": " + path);
    }
    std::string line;
    while (std::getline(in, line)) {
        on_line(line);
    }
}

GeneFeatureIndex load_gene_feature_index(
    const std::string& gtf_path,
    const std::vector<std::regex>& gene_patterns) {

    GeneFeatureIndex out;
    if (gtf_path.empty()) {
        return out;
    }

    auto keep_gene = [&](const std::string& gene_name) {
        if (gene_patterns.empty()) {
            return true;
        }
        for (const auto& pattern : gene_patterns) {
            if (std::regex_search(gene_name, pattern)) {
                return true;
            }
        }
        return false;
    };

    std::unordered_map<std::string, GeneModel> genes_by_id;

    auto ingest_line = [&](const std::string& line) {
        if (line.empty() || line[0] == '#') {
            return;
        }
        const auto fields = split_string(line, '\t');
        if (fields.size() < 9) {
            return;
        }

        const std::string feature = fields[2];
        if (feature != "gene" && feature != "exon" && feature != "CDS" && feature != "UTR") {
            return;
        }

        std::size_t start1 = 0;
        std::size_t end1 = 0;
        try {
            start1 = static_cast<std::size_t>(std::stoull(fields[3]));
            end1 = static_cast<std::size_t>(std::stoull(fields[4]));
        } catch (const std::exception&) {
            return;
        }
        if (start1 == 0 || end1 < start1) {
            return;
        }
        const std::size_t start0 = start1 - 1;
        const std::size_t end0 = end1;

        std::string gene_id = gtf_attribute_value(fields[8], "gene_id");
        std::string gene_name = gtf_attribute_value(fields[8], "gene_name");
        if (gene_name.empty()) {
            gene_name = gene_id;
        }
        if (gene_id.empty()) {
            gene_id = gene_name;
        }
        if (gene_name.empty()) {
            gene_name = ".";
        }
        if (gene_id.empty()) {
            gene_id = "GENE_UNKNOWN";
        }

        const std::string chrom = normalize_chrom(fields[0]);
        if (chrom.empty()) {
            return;
        }

        auto& model = genes_by_id[gene_id];
        if (model.chrom.empty()) {
            model.chrom = chrom;
        }
        if (model.gene_name.empty() || model.gene_name == ".") {
            model.gene_name = gene_name;
        }
        model.start0 = std::min(model.start0, start0);
        model.end0 = std::max(model.end0, end0);

        if (feature == "exon") {
            model.exons.push_back({start0, end0});
        }

        if (feature == "CDS") {
            if (keep_gene(gene_name)) {
                add_interval(out.cds, chrom, start0, end0, gene_name);
            }
        } else if (feature == "UTR") {
            if (keep_gene(gene_name)) {
                add_interval(out.utr, chrom, start0, end0, gene_name);
            }
        }
    };

    iterate_text_lines(gtf_path, ingest_line, "GTF");

    for (const auto& [gene_id, model] : genes_by_id) {
        (void)gene_id;
        if (model.chrom.empty() ||
            model.start0 == std::numeric_limits<std::size_t>::max() ||
            model.end0 <= model.start0) {
            continue;
        }
        if (!keep_gene(model.gene_name)) {
            continue;
        }

        add_interval(out.gene, model.chrom, model.start0, model.end0, model.gene_name);

        std::vector<std::pair<std::size_t, std::size_t>> clipped_exons;
        clipped_exons.reserve(model.exons.size());
        for (const auto& [s, e] : model.exons) {
            const std::size_t cs = std::max(s, model.start0);
            const std::size_t ce = std::min(e, model.end0);
            if (ce > cs) {
                clipped_exons.push_back({cs, ce});
                add_interval(out.exon, model.chrom, cs, ce, model.gene_name);
            }
        }

        const auto merged_exons = merge_intervals(std::move(clipped_exons));
        if (merged_exons.empty()) {
            continue;
        }
        std::size_t cursor = model.start0;
        for (const auto& [s, e] : merged_exons) {
            if (s > cursor) {
                add_interval(out.intron, model.chrom, cursor, s, model.gene_name);
            }
            cursor = std::max(cursor, e);
        }
        if (cursor < model.end0) {
            add_interval(out.intron, model.chrom, cursor, model.end0, model.gene_name);
        }
    }

    sort_and_unique_intervals(out.gene);
    sort_and_unique_intervals(out.exon);
    sort_and_unique_intervals(out.cds);
    sort_and_unique_intervals(out.utr);
    sort_and_unique_intervals(out.intron);

    return out;
}

struct OverlapResult {
    bool hit = false;
    std::unordered_set<std::string> genes;
};

OverlapResult query_overlap(
    const std::unordered_map<std::string, std::vector<GeneInterval>>& by_chrom,
    const std::string& raw_chrom,
    std::size_t start0,
    std::size_t end0) {

    OverlapResult out;
    if (end0 <= start0) {
        return out;
    }

    auto find_vec = [&](const std::string& key) -> const std::vector<GeneInterval>* {
        const auto it = by_chrom.find(key);
        if (it == by_chrom.end()) {
            return nullptr;
        }
        return &it->second;
    };

    const std::string chrom = normalize_chrom(raw_chrom);
    const std::vector<GeneInterval>* vec = find_vec(chrom);
    if (vec == nullptr && chrom.rfind("chr", 0) == 0 && chrom.size() > 3) {
        vec = find_vec(chrom.substr(3));
    }
    if (vec == nullptr && !chrom.empty() && chrom.rfind("chr", 0) != 0) {
        vec = find_vec("chr" + chrom);
    }
    if (vec == nullptr || vec->empty()) {
        return out;
    }

    const auto& intervals = *vec;
    auto it = std::lower_bound(
        intervals.begin(),
        intervals.end(),
        start0,
        [](const GeneInterval& interval, std::size_t value) {
            return interval.start0 < value;
        });

    std::size_t idx = static_cast<std::size_t>(std::distance(intervals.begin(), it));
    while (idx > 0 && intervals[idx - 1].end0 > start0) {
        --idx;
    }

    for (std::size_t i = idx; i < intervals.size(); ++i) {
        const auto& curr = intervals[i];
        if (curr.start0 >= end0) {
            break;
        }
        if (curr.end0 <= start0) {
            continue;
        }
        out.hit = true;
        if (!curr.gene_name.empty() && curr.gene_name != ".") {
            out.genes.insert(curr.gene_name);
        }
    }

    return out;
}

EventGeneHits annotate_event_with_genes(
    const EventRecord& ev,
    const GeneFeatureIndex* gene_index) {

    EventGeneHits out;
    if (gene_index == nullptr || gene_index->empty()) {
        return out;
    }

    if (ev.pos1 == 0) {
        return out;
    }
    const std::size_t start0 = ev.pos1 - 1;
    const std::size_t raw_end0 = (ev.end1 > ev.pos1) ? ev.end1 : ev.pos1;
    const std::size_t end0 = std::max(start0 + 1, raw_end0);

    auto absorb = [&](const OverlapResult& hit, bool& flag) {
        if (!hit.hit) {
            return;
        }
        flag = true;
        for (const auto& g : hit.genes) {
            out.genes.push_back(g);
        }
    };

    absorb(query_overlap(gene_index->gene, ev.chrom, start0, end0), out.gene);
    absorb(query_overlap(gene_index->exon, ev.chrom, start0, end0), out.exon);
    absorb(query_overlap(gene_index->cds, ev.chrom, start0, end0), out.cds);
    absorb(query_overlap(gene_index->utr, ev.chrom, start0, end0), out.utr);
    absorb(query_overlap(gene_index->intron, ev.chrom, start0, end0), out.intron);

    std::sort(out.genes.begin(), out.genes.end());
    out.genes.erase(std::unique(out.genes.begin(), out.genes.end()), out.genes.end());
    return out;
}

std::size_t size_bin_index(std::size_t length_bp, const std::vector<std::size_t>& bins) {
    const auto it = std::lower_bound(bins.begin(), bins.end(), length_bp);
    return static_cast<std::size_t>(std::distance(bins.begin(), it));
}

std::vector<std::string> size_bin_labels(const std::vector<std::size_t>& bins) {
    std::vector<std::string> labels;
    labels.reserve(bins.size() + 1);
    if (bins.empty()) {
        labels.push_back("all");
        return labels;
    }
    labels.push_back("lt_" + std::to_string(bins.front()));
    for (std::size_t i = 0; i + 1 < bins.size(); ++i) {
        labels.push_back(std::to_string(bins[i]) + "_" + std::to_string(bins[i + 1]));
    }
    labels.push_back("ge_" + std::to_string(bins.back()));
    return labels;
}

std::string feature_signature_string(const HapFeatures& f) {
    std::ostringstream out;
    out
        << f.n_events << '|'
        << (f.has_ins ? 1 : 0) << '|' << f.n_ins << '|'
        << (f.has_del ? 1 : 0) << '|' << f.n_del << '|'
        << (f.has_inv ? 1 : 0) << '|' << f.n_inv << '|'
        << (f.has_dup ? 1 : 0) << '|' << f.n_dup << '|'
        << f.n_events_gene << '|'
        << f.n_events_exon << '|'
        << f.n_events_cds << '|'
        << f.n_events_utr << '|'
        << f.n_events_intron << '|'
        << f.genes_any.size();

    auto append_vec = [&](const std::vector<std::size_t>& values) {
        out << '|';
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i > 0) {
                out << ',';
            }
            out << values[i];
        }
    };

    append_vec(f.ins_bins);
    append_vec(f.del_bins);
    append_vec(f.inv_bins);
    append_vec(f.dup_bins);

    return out.str();
}

ParsedVcf parse_region_vcf(
    const std::string& vcf_path,
    const GeneFeatureIndex* gene_index) {

    ParsedVcf out;
    std::unordered_set<std::size_t> seen_bubbles;

    auto on_line = [&](const std::string& line) {
        if (line.empty()) {
            return;
        }
        if (line.rfind("##", 0) == 0) {
            return;
        }
        if (line.rfind("#CHROM", 0) == 0) {
            const auto fields = split_string(line, '\t');
            if (fields.size() > 9) {
                out.samples.assign(fields.begin() + 9, fields.end());
            }
            return;
        }
        if (line[0] == '#') {
            return;
        }

        const auto fields = split_string(line, '\t');
        if (fields.size() < 8) {
            return;
        }

        const auto info = parse_info_field(fields[7]);
        const auto it_bubble = info.find("BUBBLE_ID");
        if (it_bubble == info.end() || it_bubble->second.empty()) {
            return;
        }

        EventRecord ev;
        ev.bubble_id = parse_size_or_throw(it_bubble->second, "BUBBLE_ID");
        ev.record_id = fields[2];
        ev.chrom = normalize_chrom(fields[0]);
        ev.pos1 = parse_size_or_throw(fields[1], "POS");
        ev.end1 = ev.pos1;
        if (const auto it_end = info.find("END"); it_end != info.end() && !it_end->second.empty()) {
            ev.end1 = parse_size_or_throw(it_end->second, "END");
        }

        if (const auto it_svtype = info.find("SVTYPE"); it_svtype != info.end()) {
            ev.svtype = it_svtype->second;
        }
        if (const auto it_event = info.find("EVENT"); it_event != info.end()) {
            ev.event_type = it_event->second;
        }
        if (const auto it_subtype = info.find("INS_SUBTYPE"); it_subtype != info.end()) {
            ev.ins_subtype = it_subtype->second;
        }
        if (ev.svtype.empty()) {
            ev.svtype = ev.event_type;
        }
        if (ev.event_type.empty()) {
            ev.event_type = ev.svtype;
        }
        if (ev.ins_subtype.empty()) {
            ev.ins_subtype = ".";
        }

        if (const auto it_svlen = info.find("SVLEN"); it_svlen != info.end()) {
            ev.svlen = parse_signed_or_default(it_svlen->second, 0);
        }

        ev.has_dup =
            (ev.ins_subtype.rfind("DUP_", 0) == 0) ||
            (info.find("DUP_REF_START") != info.end()) ||
            (info.find("DUP_SIM") != info.end());

        if (fields.size() >= 10 && !out.samples.empty()) {
            std::vector<std::string> format_keys = split_string(fields[8], ':');
            int gt_idx = -1;
            for (std::size_t i = 0; i < format_keys.size(); ++i) {
                if (format_keys[i] == "GT") {
                    gt_idx = static_cast<int>(i);
                    break;
                }
            }

            const std::size_t sample_count = std::min(out.samples.size(), fields.size() - 9);
            for (std::size_t s = 0; s < sample_count; ++s) {
                const std::string& cell = fields[9 + s];
                if (cell.empty() || cell == "." || cell == "./." || cell == ".:.") {
                    continue;
                }
                const auto values = split_string(cell, ':');
                std::string gt;
                if (gt_idx >= 0 && static_cast<std::size_t>(gt_idx) < values.size()) {
                    gt = values[static_cast<std::size_t>(gt_idx)];
                } else if (!values.empty()) {
                    gt = values[0];
                }
                if (gt_has_alt_allele(gt)) {
                    ev.carrier_indices.push_back(s);
                }
            }
        }

        ev.gene_hits = annotate_event_with_genes(ev, gene_index);

        if (seen_bubbles.insert(ev.bubble_id).second) {
            out.bubble_order.push_back(ev.bubble_id);
        }
        out.bubbles[ev.bubble_id].events.push_back(std::move(ev));
    };

    iterate_text_lines(vcf_path, on_line, "VCF");

    if (out.samples.empty()) {
        throw std::runtime_error("VCF has no sample columns: " + vcf_path);
    }

    return out;
}

void write_bubble_events_table(
    const std::filesystem::path& out_path,
    std::size_t bubble_id,
    const BubbleData& bubble,
    const std::vector<std::string>& samples) {

    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Failed to write events table: " + out_path.string());
    }

    out
        << "bubble_id\tevent_index\trecord_id\tchrom\tpos\tend\tsvtype\tevent\tins_subtype\tsvlen\thas_dup"
        << "\tcarrier_count\tcarriers"
        << "\toverlap_gene\toverlap_exon\toverlap_cds\toverlap_utr\toverlap_intron\toverlap_genes\n";

    for (std::size_t i = 0; i < bubble.events.size(); ++i) {
        const auto& ev = bubble.events[i];
        std::vector<std::string> carrier_names;
        carrier_names.reserve(ev.carrier_indices.size());
        for (const auto idx : ev.carrier_indices) {
            if (idx < samples.size()) {
                carrier_names.push_back(samples[idx]);
            }
        }

        out
            << bubble_id << '\t'
            << (i + 1) << '\t'
            << ev.record_id << '\t'
            << ev.chrom << '\t'
            << ev.pos1 << '\t'
            << ev.end1 << '\t'
            << ev.svtype << '\t'
            << ev.event_type << '\t'
            << ev.ins_subtype << '\t'
            << ev.svlen << '\t'
            << (ev.has_dup ? 1 : 0) << '\t'
            << carrier_names.size() << '\t'
            << join_csv(carrier_names) << '\t'
            << (ev.gene_hits.gene ? 1 : 0) << '\t'
            << (ev.gene_hits.exon ? 1 : 0) << '\t'
            << (ev.gene_hits.cds ? 1 : 0) << '\t'
            << (ev.gene_hits.utr ? 1 : 0) << '\t'
            << (ev.gene_hits.intron ? 1 : 0) << '\t'
            << join_csv(ev.gene_hits.genes)
            << '\n';
    }
}

void write_bubble_haplotype_table(
    const std::filesystem::path& out_path,
    std::size_t bubble_id,
    const BubbleData& bubble,
    const std::vector<std::string>& samples,
    const std::vector<std::size_t>& size_bins,
    DescribeSummary* summary_out) {

    std::unordered_map<std::string, HapFeatures> by_haplotype;
    by_haplotype.reserve(samples.size() * 2 + 1);

    const std::size_t n_bins = size_bins.size() + 1;
    for (const auto& hap : samples) {
        HapFeatures init;
        init.ins_bins.assign(n_bins, 0);
        init.del_bins.assign(n_bins, 0);
        init.inv_bins.assign(n_bins, 0);
        init.dup_bins.assign(n_bins, 0);
        by_haplotype.emplace(hap, std::move(init));
    }

    for (const auto& ev : bubble.events) {
        const std::size_t len_bp = static_cast<std::size_t>(std::llabs(ev.svlen));
        const std::size_t bin_idx = size_bin_index(len_bp, size_bins);

        for (const auto sample_idx : ev.carrier_indices) {
            if (sample_idx >= samples.size()) {
                continue;
            }
            const std::string& hap = samples[sample_idx];
            auto it = by_haplotype.find(hap);
            if (it == by_haplotype.end()) {
                continue;
            }
            auto& f = it->second;
            f.n_events += 1;

            const std::string type = ev.event_type.empty() ? ev.svtype : ev.event_type;
            if (type == "INS") {
                f.has_ins = true;
                f.n_ins += 1;
                if (bin_idx < f.ins_bins.size()) {
                    f.ins_bins[bin_idx] += 1;
                }
            } else if (type == "DEL") {
                f.has_del = true;
                f.n_del += 1;
                if (bin_idx < f.del_bins.size()) {
                    f.del_bins[bin_idx] += 1;
                }
            } else if (type == "INV") {
                f.has_inv = true;
                f.n_inv += 1;
                if (bin_idx < f.inv_bins.size()) {
                    f.inv_bins[bin_idx] += 1;
                }
            }

            if (ev.has_dup) {
                f.has_dup = true;
                f.n_dup += 1;
                if (bin_idx < f.dup_bins.size()) {
                    f.dup_bins[bin_idx] += 1;
                }
            }

            if (ev.gene_hits.gene) {
                f.n_events_gene += 1;
            }
            if (ev.gene_hits.exon) {
                f.n_events_exon += 1;
            }
            if (ev.gene_hits.cds) {
                f.n_events_cds += 1;
            }
            if (ev.gene_hits.utr) {
                f.n_events_utr += 1;
            }
            if (ev.gene_hits.intron) {
                f.n_events_intron += 1;
            }
            for (const auto& gene_name : ev.gene_hits.genes) {
                f.genes_any.insert(gene_name);
            }
        }
    }

    std::vector<std::string> haps_sorted = samples;
    std::sort(haps_sorted.begin(), haps_sorted.end());

    std::unordered_map<std::string, std::size_t> signature_counts;
    signature_counts.reserve(haps_sorted.size() * 2 + 1);
    for (const auto& hap : haps_sorted) {
        const auto it = by_haplotype.find(hap);
        if (it == by_haplotype.end()) {
            continue;
        }
        const std::string key = feature_signature_string(it->second);
        signature_counts[key] += 1;
    }

    std::vector<std::pair<std::string, std::size_t>> signature_vec(signature_counts.begin(), signature_counts.end());
    std::sort(signature_vec.begin(), signature_vec.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) {
            return a.second > b.second;
        }
        return a.first < b.first;
    });

    std::unordered_map<std::string, std::string> signature_ids;
    signature_ids.reserve(signature_vec.size() * 2 + 1);
    for (std::size_t i = 0; i < signature_vec.size(); ++i) {
        signature_ids[signature_vec[i].first] = "SIG" + std::to_string(i + 1);
    }

    const auto labels = size_bin_labels(size_bins);

    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Failed to write haplotype feature table: " + out_path.string());
    }

    out
        << "bubble_id\thaplotype"
        << "\tn_events\thas_ins\tn_ins\thas_del\tn_del\thas_inv\tn_inv\thas_dup\tn_dup"
        << "\tn_events_gene\tn_events_exon\tn_events_cds\tn_events_utr\tn_events_intron\tn_genes_any";

    auto emit_bin_headers = [&](const std::string& prefix) {
        for (const auto& label : labels) {
            out << '\t' << "n_" << prefix << "_" << label;
        }
    };
    emit_bin_headers("ins_len");
    emit_bin_headers("del_len");
    emit_bin_headers("inv_len");
    emit_bin_headers("dup_len");

    out << "\tfeature_signature\tsignature_group_size\n";

    auto emit_bins = [&](const std::vector<std::size_t>& values) {
        for (const auto v : values) {
            out << '\t' << v;
        }
    };

    for (const auto& hap : haps_sorted) {
        const auto it = by_haplotype.find(hap);
        if (it == by_haplotype.end()) {
            continue;
        }
        const auto& f = it->second;
        const std::string sig_key = feature_signature_string(f);
        const auto it_sig_id = signature_ids.find(sig_key);
        const auto it_sig_count = signature_counts.find(sig_key);
        const std::string sig_id = (it_sig_id == signature_ids.end()) ? "SIG0" : it_sig_id->second;
        const std::size_t sig_count = (it_sig_count == signature_counts.end()) ? 0 : it_sig_count->second;

        out
            << bubble_id << '\t'
            << hap << '\t'
            << f.n_events << '\t'
            << (f.has_ins ? 1 : 0) << '\t'
            << f.n_ins << '\t'
            << (f.has_del ? 1 : 0) << '\t'
            << f.n_del << '\t'
            << (f.has_inv ? 1 : 0) << '\t'
            << f.n_inv << '\t'
            << (f.has_dup ? 1 : 0) << '\t'
            << f.n_dup << '\t'
            << f.n_events_gene << '\t'
            << f.n_events_exon << '\t'
            << f.n_events_cds << '\t'
            << f.n_events_utr << '\t'
            << f.n_events_intron << '\t'
            << f.genes_any.size();

        emit_bins(f.ins_bins);
        emit_bins(f.del_bins);
        emit_bins(f.inv_bins);
        emit_bins(f.dup_bins);

        out
            << '\t' << sig_id
            << '\t' << sig_count
            << '\n';

        if (summary_out != nullptr) {
            summary_out->haplotype_rows += 1;
        }
    }
}

std::vector<std::regex> compile_gene_regexes(const std::vector<std::string>& patterns) {
    std::vector<std::regex> compiled;
    compiled.reserve(patterns.size());
    for (const auto& raw : patterns) {
        if (raw.empty()) {
            continue;
        }
        try {
            compiled.emplace_back(raw, std::regex::icase);
        } catch (const std::regex_error& e) {
            throw std::runtime_error("Invalid regex in --gene-match: '" + raw + "' (" + e.what() + ")");
        }
    }
    return compiled;
}

} // namespace

void describe_from_region_vcf(
    const DescribeOptions& options,
    DescribeSummary* summary_out) {

    if (options.vcf_in_path.empty()) {
        throw std::runtime_error("describe requires --vcf-in");
    }

    DescribeSummary local_summary;

    std::filesystem::create_directories(options.out_dir);

    const auto gene_patterns = compile_gene_regexes(options.gene_match_patterns);
    const GeneFeatureIndex gene_index = load_gene_feature_index(options.gtf_path, gene_patterns);
    const GeneFeatureIndex* gene_ptr = gene_index.empty() ? nullptr : &gene_index;

    const ParsedVcf parsed = parse_region_vcf(options.vcf_in_path, gene_ptr);

    std::vector<std::size_t> bins = options.size_bins;
    if (bins.empty()) {
        bins = {100, 1000};
    }
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    if (bins.empty() || bins.front() == 0) {
        throw std::runtime_error("size bins must contain positive integers");
    }

    const std::filesystem::path index_path = std::filesystem::path(options.out_dir) / "describe.index.tsv";
    std::ofstream index_out(index_path);
    if (!index_out) {
        throw std::runtime_error("Failed to write describe index: " + index_path.string());
    }
    index_out << "bubble_id\tevents_table\thaplotype_table\tn_events\tn_haplotypes\n";

    for (const auto bubble_id : parsed.bubble_order) {
        const auto it = parsed.bubbles.find(bubble_id);
        if (it == parsed.bubbles.end()) {
            continue;
        }
        const auto& bubble = it->second;

        const std::filesystem::path events_path =
            std::filesystem::path(options.out_dir) / ("bubble_" + std::to_string(bubble_id) + ".events.tsv");
        const std::filesystem::path hap_path =
            std::filesystem::path(options.out_dir) / ("bubble_" + std::to_string(bubble_id) + ".haplotype_features.tsv");

        write_bubble_events_table(events_path, bubble_id, bubble, parsed.samples);

        const std::size_t rows_before = local_summary.haplotype_rows;
        write_bubble_haplotype_table(hap_path, bubble_id, bubble, parsed.samples, bins, &local_summary);
        const std::size_t rows_after = local_summary.haplotype_rows;

        index_out
            << bubble_id << '\t'
            << events_path.filename().string() << '\t'
            << hap_path.filename().string() << '\t'
            << bubble.events.size() << '\t'
            << (rows_after - rows_before)
            << '\n';

        local_summary.bubbles += 1;
        local_summary.events += bubble.events.size();
        local_summary.files_written += 2;

        if (!options.quiet) {
            std::cerr
                << "[describe] bubble " << bubble_id
                << " events=" << bubble.events.size()
                << " hap_rows=" << (rows_after - rows_before)
                << "\n";
        }
    }

    local_summary.files_written += 1; // index table

    if (summary_out != nullptr) {
        *summary_out = local_summary;
    }
}

} // namespace panvar
