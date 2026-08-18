#include "panvar/variant_nodes.hpp"

#include "panvar/bubbles.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace panvar {

// `bubble_members` is every node each bubble owns (interior plus both boundaries). A call's nodes are
// the whole basis on which a truth event is attributed to it, so a node that is not in the graph, or
// belongs to a different site, silently turns a called event into a missed one -- the input error and
// the caller failure produce the identical output. Nothing here is skipped or repaired: an input that
// cannot be trusted is refused.
VariantNodes load_variant_nodes(
        const std::string& path,
        const std::unordered_map<std::size_t, std::unordered_set<std::string>>& bubble_members,
        const std::string& module) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error(module + ": failed to open variant nodes: " + path);
    auto split_tab = [](const std::string& s) {
        std::vector<std::string> f;
        std::string cur;
        for (char c : s) { if (c == '\t') { f.push_back(cur); cur.clear(); } else cur += c; }
        f.push_back(cur);
        return f;
    };
    VariantNodes cv;
    std::string line;
    std::size_t ncol = 0, lineno = 0;
    std::unordered_set<std::string> seen_ids;
    static const char* kHead[] = {"variant_id", "bubble_id", "svtype", "node_ids"};
    while (std::getline(in, line)) {
        ++lineno;
        if (line.empty()) continue;
        const std::vector<std::string> f = split_tab(line);
        if (ncol == 0) {
            // A missing header is not cosmetic: the first data row would be eaten as one, dropping a
            // real call and reporting its event as missed.
            if (f.size() < 4)
                throw std::runtime_error(module + ": " + path + ": expected a header with at least 4 tab-separated "
                                         "columns, found " + std::to_string(f.size()));
            for (int i = 0; i < 4; ++i)
                if (f[static_cast<std::size_t>(i)] != kHead[i])
                    throw std::runtime_error(module + ": " + path + ": column " + std::to_string(i + 1) + " of the header is '" +
                                             f[static_cast<std::size_t>(i)] + "', expected '" + kHead[i] +
                                             "'. This is not a variant_nodes.tsv written by `call`");
            ncol = f.size();
            continue;
        }
        if (f.size() != ncol)
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": has " +
                                     std::to_string(f.size()) + " fields, the header has " +
                                     std::to_string(ncol));
        if (f[0].empty())
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": empty variant_id");
        if (!seen_ids.insert(f[0]).second)
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": duplicate variant_id '" +
                                     f[0] + "'; which record a truth event is attributed to would be undefined");
        std::size_t bid = 0;
        try {
            std::size_t used = 0;
            bid = static_cast<std::size_t>(std::stoull(f[1], &used));
            if (used != f[1].size()) throw std::invalid_argument("trailing");
        } catch (const std::exception&) {
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": bubble_id '" + f[1] +
                                     "' is not a number");
        }
        const auto mit = bubble_members.find(bid);
        if (mit == bubble_members.end())
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": bubble id " +
                                     std::to_string(bid) + " is not in the bubbles CSV; these outputs "
                                     "are from different runs");
        CalledRecord rec;
        rec.variant_id = f[0];
        rec.svtype = f[2];
        auto& un = cv.nodes[bid];
        std::string tok;
        auto take = [&]() {
            if (tok.empty()) return;
            if (!mit->second.count(tok))
                throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": record '" + f[0] +
                                         "' names node '" + tok + "', which bubble " + std::to_string(bid) +
                                         " does not contain. A stale node turns a called event into a "
                                         "missed one with nothing to show for it");
            rec.nodes.insert(tok);
            un.insert(tok);
            tok.clear();
        };
        for (char c : f[3]) { if (c == ',') take(); else tok += c; }
        take();
        if (rec.nodes.empty())
            throw std::runtime_error(module + ": " + path + ":" + std::to_string(lineno) + ": record '" + f[0] +
                                     "' names no nodes");
        cv.bubble_ids.insert(bid);
        cv.svtypes[bid].insert(f[2]);
        cv.records[bid].push_back(std::move(rec));
    }
    if (ncol == 0) throw std::runtime_error(module + ": " + path + ": file is empty");
    // Deterministic record order, so which record a truth event is attributed to cannot depend on the
    // order rows happen to sit in the file.
    for (auto& [bid, recs] : cv.records)
        std::sort(recs.begin(), recs.end(),
                  [](const CalledRecord& a, const CalledRecord& b) { return a.variant_id < b.variant_id; });
    return cv;
}


std::unordered_map<std::size_t, std::unordered_set<std::string>> bubble_member_nodes(
        const std::vector<Bubble>& bubbles) {
    std::unordered_map<std::size_t, std::unordered_set<std::string>> out;
    for (const Bubble& b : bubbles) {
        std::unordered_set<std::string> members;
        members.insert(b.source);
        members.insert(b.sink);
        for (const std::string& n : b.inside) members.insert(n);
        out.emplace(b.id, std::move(members));
    }
    return out;
}

} // namespace panvar
