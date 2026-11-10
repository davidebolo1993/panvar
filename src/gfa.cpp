#include "panvar/gfa.hpp"

#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#include <utility>

namespace panvar {
namespace {

struct RawEdge {
    std::string from;
    char from_orient = '+';
    std::string to;
    char to_orient = '+';
    int overlap = 0;
    std::size_t line_no = 0;
};

std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> parts;
    std::string item;
    std::stringstream ss(s);
    while (std::getline(ss, item, delim)) {
        parts.push_back(item);
    }
    return parts;
}

bool neighbor_exists(const std::vector<Neighbor>& vec, const Neighbor& n) {
    return std::any_of(vec.begin(), vec.end(), [&](const Neighbor& x) {
        return x.node_id == n.node_id && x.side == n.side && x.overlap == n.overlap;
    });
}

void add_neighbor(std::vector<Neighbor>& vec, const Neighbor& n) {
    if (!neighbor_exists(vec, n)) {
        vec.push_back(n);
    }
}

std::vector<PathStep> parse_p_steps(const std::string& field) {
    std::vector<PathStep> steps;
    // getline drops a trailing empty field, so "1+,2+," would lose its empty last step silently.
    if (!field.empty() && field.back() == ',') {
        throw std::runtime_error("Invalid P-line step list: it ends with a comma, an empty step");
    }
    for (const auto& token : split(field, ',')) {
        // An empty step ("1+,,2+") or one with no id ("+") is refused: skipping it would parse a damaged
        // path one step short.
        if (token.size() < 2) {
            throw std::runtime_error("Invalid P-line step '" + token + "': empty, or no segment id");
        }
        const char orient = token.back();
        if (orient != '+' && orient != '-') {
            throw std::runtime_error("Invalid P-line segment orientation: " + token);
        }
        PathStep step;
        step.node_id = token.substr(0, token.size() - 1);
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

// `*` on an L line means the overlap is UNKNOWN, not that it is zero. Reporting it as 0 lets a caller
// that concatenates segments believe it verified something it did not, so it is distinguished here and
// consumers decide: -1 = unknown.
// Only '*' and a plain match length ("<n>M") are understood; anything else ("0M2I", "0", "-1M") is
// refused rather than read by its leading digits.
int parse_overlap(const std::string& field, std::size_t line_no) {
    if (field == "*") {
        return -1;
    }
    const bool digits = field.size() >= 2 && field.size() <= 10 && field.back() == 'M' &&
                        std::all_of(field.begin(), field.end() - 1,
                                    [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
    if (!digits) {
        throw std::runtime_error("Invalid L line at " + std::to_string(line_no) + ": overlap '" + field +
                                 "' is neither '*' nor a match length like 0M");
    }
    return std::stoi(field.substr(0, field.size() - 1));
}

} // namespace

std::vector<std::string> Graph::neighbors_of(const std::string& node_id) const {
    const auto it = nodes.find(node_id);
    if (it == nodes.end()) {
        return {};
    }
    std::unordered_set<std::string> uniq;
    for (const auto& n : it->second.start) {
        uniq.insert(n.node_id);
    }
    for (const auto& n : it->second.end) {
        uniq.insert(n.node_id);
    }
    std::vector<std::string> out(uniq.begin(), uniq.end());
    std::sort(out.begin(), out.end());
    return out;
}

bool Graph::has_edge_between(const std::string& a, const std::string& b) const {
    const auto it = nodes.find(a);
    if (it == nodes.end()) {
        return false;
    }
    auto contains = [&](const std::vector<Neighbor>& v) {
        return std::any_of(v.begin(), v.end(), [&](const Neighbor& n) {
            return n.node_id == b;
        });
    };
    return contains(it->second.start) || contains(it->second.end);
}

void check_p_overlaps(const std::string& field, std::size_t line_no) {
    if (field.empty() || field == "*") return;
    for (const auto& token : split(field, ',')) {
        if (token == "*") continue;
        const bool zero = token.size() >= 2 && token.back() == 'M' &&
                          std::all_of(token.begin(), token.end() - 1, [](char c) { return c == '0'; });
        if (!zero)
            throw std::runtime_error("Invalid P line at " + std::to_string(line_no) + ": overlap '" + token +
                                     "' is not 0M or '*' (paths are spelled by concatenating their segments)");
    }
}

std::string w_literal_as_path_name(const std::string& query) {
    const std::size_t colon = query.rfind(':');
    const std::size_t dash = query.rfind('-');
    if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1 || dash + 1 >= query.size())
        return std::string();
    const std::string a = query.substr(colon + 1, dash - colon - 1), b = query.substr(dash + 1);
    const auto digits = [](const std::string& t) {
        return !t.empty() && t.size() < 19 &&
               std::all_of(t.begin(), t.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
    };
    if (!digits(a) || (!digits(b) && b != "*")) return std::string();
    return query.substr(0, colon + 1) + std::to_string(std::stoull(a) + 1) + "-" + b;
}

std::string w_path_name(const std::string& sample, const std::string& hap, const std::string& seqid,
                        const std::string& start, const std::string& end) {
    const auto numeric = [](const std::string& t) {
        return !t.empty() && t.size() < 19 &&
               std::all_of(t.begin(), t.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
    };
    // No start ('*', as GFA 1.1 allows): no interval in the name, as "seqid:*-*" would read as a contig of
    // that name. A start without an end keeps the start, written "start-*", since positions need only the
    // start.
    if (!numeric(start)) return sample + "#" + hap + "#" + seqid;
    return sample + "#" + hap + "#" + seqid + ":" + std::to_string(std::stoull(start) + 1) + "-" +
           (numeric(end) ? end : std::string("*"));
}

Graph parse_gfa(const std::string& gfa_path, const ParseGfaOptions& options) {
    GzLineReader in(gfa_path);
    if (!in.ok()) {
        throw std::runtime_error("Failed to open GFA file: " + gfa_path);
    }

    Graph graph;
    std::vector<RawEdge> raw_edges;

    std::string line;
    std::size_t line_no = 0;
    while (in.getline(line)) {
        ++line_no;
        if (line.empty()) {
            continue;
        }
        const auto fields = split(line, '\t');
        if (fields.empty()) {
            continue;
        }

        const std::string& rec = fields[0];
        if (rec == "S") {
            if (fields.size() < 3) {
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no));
            }
            if (fields[1].empty()) {
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no) + ": empty segment id");
            }
            // A second S line for one id is refused: keeping either sequence would leave every link and path
            // valid, and two GFAs concatenated with numeric ids would read as one chimeric graph.
            if (graph.nodes.count(fields[1])) {
                throw std::runtime_error("Invalid S line at " + std::to_string(line_no) + ": segment '" +
                                         fields[1] + "' is defined more than once");
            }
            Node& node = graph.nodes[fields[1]];
            node.id = fields[1];
            if (options.include_sequences && fields[2] != "*") {
                node.sequence = fields[2];
            } else {
                node.sequence.clear();
            }
        } else if (rec == "L") {
            if (fields.size() < 6) {
                throw std::runtime_error("Invalid L line at " + std::to_string(line_no));
            }
            RawEdge edge;
            edge.from = fields[1];
            if ((fields[2] != "+" && fields[2] != "-") || (fields[4] != "+" && fields[4] != "-")) {
                throw std::runtime_error("Invalid L line at " + std::to_string(line_no) +
                                         ": orientations must be '+' or '-'");
            }
            edge.from_orient = fields[2].empty() ? '+' : fields[2][0];
            edge.to = fields[3];
            edge.to_orient = fields[4].empty() ? '+' : fields[4][0];
            edge.overlap = parse_overlap(fields[5], line_no);
            edge.line_no = line_no;
            raw_edges.push_back(std::move(edge));
        } else if (rec == "P" && options.include_paths) {
            if (fields.size() < 3) {
                throw std::runtime_error("Invalid P line at " + std::to_string(line_no));
            }
            PathRecord p;
            p.name = fields[1];
            p.type = 'P';
            p.steps = parse_p_steps(fields[2]);
            if (fields.size() >= 4) check_p_overlaps(fields[3], line_no);
            graph.paths.push_back(std::move(p));
        } else if (rec == "W" && options.include_paths) {
            if (fields.size() < 7) {
                throw std::runtime_error("Invalid W line at " + std::to_string(line_no));
            }
            PathRecord w;
            w.name = w_path_name(fields[1], fields[2], fields[3], fields[4], fields[5]);
            w.type = 'W';
            w.steps = parse_w_steps(fields[6], line_no);
            graph.paths.push_back(std::move(w));
        }
    }

    for (const auto& edge : raw_edges) {
        // A link to a segment the file does not define is damage, refused rather than dropped unseen.
        for (const std::string* id : {&edge.from, &edge.to})
            if (graph.nodes.find(*id) == graph.nodes.end())
                throw std::runtime_error("Invalid L line at " + std::to_string(edge.line_no) + ": segment " + *id +
                                         " is not in the graph");

        const bool from_start = (edge.from_orient == '-');
        const bool to_end = (edge.to_orient == '-');

        if (from_start && to_end) {
            add_neighbor(graph.nodes[edge.from].start, Neighbor{edge.to, 1, edge.overlap});
            add_neighbor(graph.nodes[edge.to].end, Neighbor{edge.from, 0, edge.overlap});
        } else if (from_start && !to_end) {
            add_neighbor(graph.nodes[edge.from].start, Neighbor{edge.to, 0, edge.overlap});
            add_neighbor(graph.nodes[edge.to].start, Neighbor{edge.from, 0, edge.overlap});
        } else if (!from_start && !to_end) {
            add_neighbor(graph.nodes[edge.from].end, Neighbor{edge.to, 0, edge.overlap});
            add_neighbor(graph.nodes[edge.to].start, Neighbor{edge.from, 1, edge.overlap});
        } else {
            add_neighbor(graph.nodes[edge.from].end, Neighbor{edge.to, 1, edge.overlap});
            add_neighbor(graph.nodes[edge.to].end, Neighbor{edge.from, 1, edge.overlap});
        }
    }

    return graph;
}

} // namespace panvar
