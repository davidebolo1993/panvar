#include "panvar/gfa_io.hpp"

#include "panvar/gz_reader.hpp"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace panvar {
namespace {

std::vector<std::string> split_tab(const std::string& s) {
    std::vector<std::string> parts;
    std::string item;
    std::stringstream ss(s);
    while (std::getline(ss, item, '\t')) {
        parts.push_back(item);
    }
    return parts;
}

std::vector<PathStep> parse_p_steps(const std::string& field) {
    std::vector<PathStep> steps;
    std::string item;
    std::stringstream ss(field);
    while (std::getline(ss, item, ',')) {
        if (item.empty()) {
            continue;
        }
        const char orient = item.back();
        if (orient != '+' && orient != '-') {
            throw std::runtime_error("Invalid P-line segment orientation: " + item);
        }
        PathStep step;
        step.node_id = item.substr(0, item.size() - 1);
        step.reverse = (orient == '-');
        steps.push_back(std::move(step));
    }
    return steps;
}

std::vector<PathStep> parse_w_steps(const std::string& walk) {
    std::vector<PathStep> steps;
    std::size_t i = 0;
    while (i < walk.size()) {
        if (walk[i] != '>' && walk[i] != '<') {
            ++i;
            continue;
        }
        const bool reverse = (walk[i] == '<');
        ++i;
        const std::size_t start = i;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') {
            ++i;
        }
        if (i == start) {
            continue;
        }
        PathStep step;
        step.node_id = walk.substr(start, i - start);
        step.reverse = reverse;
        steps.push_back(std::move(step));
    }
    return steps;
}

std::string format_p_steps(const std::vector<PathStep>& steps) {
    std::string out;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (i > 0) {
            out.push_back(',');
        }
        out += steps[i].node_id;
        out.push_back(steps[i].reverse ? '-' : '+');
    }
    return out;
}

std::string format_w_steps(const std::vector<PathStep>& steps) {
    std::string out;
    out.reserve(steps.size() * 4);
    for (const PathStep& s : steps) {
        out.push_back(s.reverse ? '<' : '>');
        out += s.node_id;
    }
    return out;
}

} // namespace

GfaModel read_gfa_model(const std::string& gfa_path) {
    GzLineReader in(gfa_path);
    if (!in.ok()) {
        throw std::runtime_error("Failed to open GFA file: " + gfa_path);
    }

    GfaModel model;
    std::string line;
    std::size_t line_no = 0;
    while (in.getline(line)) {
        ++line_no;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        const char rec = line[0];
        if (rec == 'H') {
            model.header_lines.push_back(line);
            continue;
        }
        const auto fields = split_tab(line);
        if (fields.empty()) {
            continue;
        }
        if (fields[0] == "S") {
            if (fields.size() < 3) {
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no));
            }
            const std::string& id = fields[1];
            if (!model.has_node(id)) {
                model.node_order.push_back(id);
            }
            model.seq[id] = (fields[2] == "*") ? std::string() : fields[2];
        } else if (fields[0] == "L") {
            if (fields.size() < 6) {
                throw std::runtime_error("Invalid L line at " + std::to_string(line_no));
            }
            GfaEdge e;
            e.from = fields[1];
            e.from_orient = fields[2].empty() ? '+' : fields[2][0];
            e.to = fields[3];
            e.to_orient = fields[4].empty() ? '+' : fields[4][0];
            e.overlap = fields[5];
            model.edges.push_back(std::move(e));
        } else if (fields[0] == "P") {
            if (fields.size() < 3) {
                throw std::runtime_error("Invalid P line at " + std::to_string(line_no));
            }
            GfaPath p;
            p.type = 'P';
            p.name = fields[1];
            p.steps = parse_p_steps(fields[2]);
            p.overlaps = (fields.size() >= 4 && !fields[3].empty()) ? fields[3] : "*";
            model.paths.push_back(std::move(p));
        } else if (fields[0] == "W") {
            if (fields.size() < 7) {
                throw std::runtime_error("Invalid W line at " + std::to_string(line_no));
            }
            GfaPath w;
            w.type = 'W';
            w.sample = fields[1];
            w.hap = fields[2];
            w.seqid = fields[3];
            w.start = fields[4];
            w.end = fields[5];
            w.steps = parse_w_steps(fields[6]);
            model.paths.push_back(std::move(w));
        }
    }
    return model;
}

void write_gfa_model(const std::string& out_path, const GfaModel& model) {
    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Failed to write GFA file: " + out_path);
    }
    if (model.header_lines.empty()) {
        out << "H\tVN:Z:1.1\n";
    } else {
        for (const std::string& h : model.header_lines) {
            out << h << '\n';
        }
    }
    for (const std::string& id : model.node_order) {
        const auto it = model.seq.find(id);
        const std::string& s = it != model.seq.end() ? it->second : std::string();
        out << "S\t" << id << '\t' << (s.empty() ? "*" : s) << '\n';
    }
    for (const GfaEdge& e : model.edges) {
        out << "L\t" << e.from << '\t' << e.from_orient << '\t'
            << e.to << '\t' << e.to_orient << '\t' << e.overlap << '\n';
    }
    for (const GfaPath& p : model.paths) {
        if (p.type == 'W') {
            out << "W\t" << p.sample << '\t' << p.hap << '\t' << p.seqid << '\t'
                << p.start << '\t' << p.end << '\t' << format_w_steps(p.steps) << '\n';
        } else {
            out << "P\t" << p.name << '\t' << format_p_steps(p.steps) << '\t'
                << (p.overlaps.empty() ? "*" : p.overlaps) << '\n';
        }
    }
}

std::string add_new_node(GfaModel& model, const std::string& sequence) {
    std::uint64_t max_id = 0;
    for (const std::string& id : model.node_order) {
        // ids are numeric in these graphs; fall back to size-based ids otherwise.
        char* endp = nullptr;
        const std::uint64_t v = std::strtoull(id.c_str(), &endp, 10);
        if (endp != id.c_str() && *endp == '\0' && v > max_id) {
            max_id = v;
        }
    }
    std::string new_id = std::to_string(max_id + 1);
    while (model.has_node(new_id)) {
        max_id += 1;
        new_id = std::to_string(max_id + 1);
    }
    model.node_order.push_back(new_id);
    model.seq[new_id] = sequence;
    return new_id;
}

} // namespace panvar
