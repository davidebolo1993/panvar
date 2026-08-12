#include "panvar/gfa.hpp"

#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace panvar {
namespace {

struct RawEdge {
    std::string from;
    char from_orient = '+';
    std::string to;
    char to_orient = '+';
    int overlap = 0;
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
    for (const auto& token : split(field, ',')) {
        if (token.empty()) {
            continue;
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

// `*` on an L line means the overlap is UNKNOWN, not that it is zero. Reporting it as 0 lets a caller
// that concatenates segments believe it verified something it did not, so it is distinguished here and
// consumers decide: -1 = unknown.
int parse_overlap(const std::string& field) {
    if (field == "*") {
        return -1;
    }
    if (!field.empty() && field.back() == 'M') {
        return std::stoi(field.substr(0, field.size() - 1));
    }
    return std::stoi(field);
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
            edge.from_orient = fields[2].empty() ? '+' : fields[2][0];
            edge.to = fields[3];
            edge.to_orient = fields[4].empty() ? '+' : fields[4][0];
            edge.overlap = parse_overlap(fields[5]);
            raw_edges.push_back(std::move(edge));
        } else if (rec == "P" && options.include_paths) {
            if (fields.size() < 3) {
                throw std::runtime_error("Invalid P line at " + std::to_string(line_no));
            }
            PathRecord p;
            p.name = fields[1];
            p.type = 'P';
            p.steps = parse_p_steps(fields[2]);
            graph.paths.push_back(std::move(p));
        } else if (rec == "W" && options.include_paths) {
            if (fields.size() < 7) {
                throw std::runtime_error("Invalid W line at " + std::to_string(line_no));
            }
            PathRecord w;
            w.name = fields[1] + "#" + fields[2] + "#" + fields[3] + ":" + fields[4] + "-" + fields[5];
            w.type = 'W';
            w.steps = parse_w_steps(fields[6]);
            graph.paths.push_back(std::move(w));
        }
    }

    for (const auto& edge : raw_edges) {
        if (graph.nodes.find(edge.from) == graph.nodes.end() ||
            graph.nodes.find(edge.to) == graph.nodes.end()) {
            continue;
        }

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
