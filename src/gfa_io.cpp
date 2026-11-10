#include "panvar/gfa_io.hpp"

#include "panvar/cli_utils.hpp"

#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
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
    // getline drops an empty last item, so a trailing comma would vanish unread
    if (!field.empty() && field.back() == ',')
        throw std::runtime_error("Invalid P-line step list: it ends with a comma, an empty step");
    std::vector<PathStep> steps;
    std::string item;
    std::stringstream ss(field);
    while (std::getline(ss, item, ',')) {
        if (item.size() < 2) {   // as in parse_gfa: an empty step is refused, not skipped
            throw std::runtime_error("Invalid P-line step '" + item + "': empty, or no segment id");
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

// Every step is '>' or '<' then a segment id. Text before the first marker, or a marker with no id, is
// refused: skipping it would parse a damaged walk ("9>1>3") one step short.
std::vector<PathStep> parse_w_steps(const std::string& walk, std::size_t line_no) {
    std::vector<PathStep> steps;
    std::size_t i = 0;
    while (i < walk.size()) {
        if (walk[i] != '>' && walk[i] != '<') {
            throw std::runtime_error("Invalid W line at " + std::to_string(line_no) +
                                     ": walk text '" + walk.substr(i, 20) + "' is not a step");
        }
        const bool reverse = (walk[i] == '<');
        ++i;
        const std::size_t start = i;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') {
            ++i;
        }
        if (i == start) {
            throw std::runtime_error("Invalid W line at " + std::to_string(line_no) +
                                     ": a walk step has no segment id");
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
            if (id.empty()) {
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no) + ": empty segment id");
            }
            if (model.has_node(id)) {   // refused, as in parse_gfa
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no) + ": segment '" + id +
                                         "' is defined more than once");
            }
            model.node_order.push_back(id);
            model.seq[id] = (fields[2] == "*") ? std::string() : fields[2];
        } else if (fields[0] == "L") {
            if (fields.size() < 6) {
                throw std::runtime_error("Invalid L line at " + std::to_string(line_no));
            }
            GfaEdge e;
            e.from = fields[1];
            if ((fields[2] != "+" && fields[2] != "-") || (fields[4] != "+" && fields[4] != "-")) {
                throw std::runtime_error("Invalid L line at " + std::to_string(line_no) +
                                         ": orientations must be '+' or '-'");
            }
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
            if (fields.size() >= 4) check_p_overlaps(fields[3], line_no);
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
            w.steps = parse_w_steps(fields[6], line_no);
            model.paths.push_back(std::move(w));
        }
    }
    return model;
}

void write_gfa_model(const std::string& out_path, const GfaModel& model) {
    errno = 0;
    std::ofstream out(out_path);
    if (!out) {
        // close_output names the destination, not a staging file; its hint is for a failed write, so a
        // failed open gives the system's reason instead (permission denied, no such directory).
        const int err = errno;
        try {
            cli::close_output(out, out_path);
        } catch (const std::runtime_error& e) {
            std::string msg = e.what();
            const std::size_t hint = msg.find(" (disk full or I/O error?)");
            if (err != 0 && hint != std::string::npos)
                msg.replace(hint, std::string::npos, std::string(": ") + std::strerror(err));
            throw std::runtime_error(msg);
        }
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
            // A P line's overlap list has one entry per link. A module that rewrites the steps (panphorte
            // folding, refine) leaves the input's list behind, and a 4-step walk written with 5 overlaps
            // is invalid GFA; such a list is written as `*` (unspecified).
            std::string overlaps = p.overlaps.empty() ? "*" : p.overlaps;
            if (overlaps != "*") {
                const std::size_t entries = static_cast<std::size_t>(std::count(overlaps.begin(), overlaps.end(), ',')) + 1;
                if (p.steps.empty() || entries != p.steps.size() - 1) overlaps = "*";
            }
            out << "P\t" << p.name << '\t' << format_p_steps(p.steps) << '\t' << overlaps << '\n';
        }
    }
    cli::close_output(out, out_path);
}

std::string add_new_node(GfaModel& model, const std::string& sequence) {
    // A node_order that shrank was rebuilt (a sort renumbers it): scan it whole again.
    if (model.max_id_scanned > model.node_order.size()) {
        model.max_id_scanned = 0;
        model.max_id_seen = 0;
    }
    std::uint64_t max_id = model.max_id_seen;
    for (std::size_t i = model.max_id_scanned; i < model.node_order.size(); ++i) {
        const std::string& id = model.node_order[i];
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
    model.max_id_seen = std::max<std::uint64_t>(max_id, max_id + 1);
    model.max_id_scanned = model.node_order.size();
    return new_id;
}

} // namespace panvar
