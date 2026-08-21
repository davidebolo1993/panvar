#pragma once

// Reader for `call`'s <prefix>.variant_nodes.tsv, shared by every consumer.
//
// The file is one row per emitted RECORD holding the union of nodes over every merged event and every
// carrier. It has no haplotype column: two carriers of one merged record traverse different node sets,
// so which haplotype carries what is only knowable from the VCF.
//
// It is read strictly, because a call's node set is the whole basis on which anything downstream is
// attributed to that call. A stale node, a dropped row or a mis-parsed bubble id would not fail --
// it produces exactly the output a genuine caller miss produces, and the two are indistinguishable
// afterwards. So nothing here is skipped or repaired.

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {

// One emitted call. Kept whole: folding a bubble's calls into a single node set cannot distinguish
// "this event matched THAT record" from "this event touches something some record also touches".
struct CalledRecord {
    std::string variant_id;
    std::string svtype;
    std::unordered_set<std::string> nodes;
};

struct VariantNodes {
    std::unordered_set<std::size_t> bubble_ids;                             // bubbles with >=1 call
    std::unordered_map<std::size_t, std::vector<CalledRecord>> records;     // bubble -> records
    std::unordered_map<std::size_t, std::unordered_set<std::string>> nodes; // bubble -> node union
    std::map<std::size_t, std::set<std::string>> svtypes;                   // bubble -> svtypes
};

// `bubble_members` is every node each bubble owns (interior plus both boundaries), keyed by bubble id.
// Every row is checked against it: the bubble must be present, and every node the row names must
// belong to that bubble. `module` prefixes error messages.
//
// Refuses: a header that is not `call`'s, a row whose field count differs from the header's, an empty
// or duplicate variant_id, a non-numeric bubble id, a bubble id absent from `bubble_members`, a node
// the named bubble does not contain, and a row naming no nodes.
VariantNodes load_variant_nodes(
    const std::string& path,
    const std::unordered_map<std::size_t, std::unordered_set<std::string>>& bubble_members,
    const std::string& module);

// Every node each bubble owns, for the map above.
struct Bubble;
std::unordered_map<std::size_t, std::unordered_set<std::string>> bubble_member_nodes(
    const std::vector<Bubble>& bubbles);

} // namespace panvar
