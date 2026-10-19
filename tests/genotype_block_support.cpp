// M3-E gates: the shared full-data-likelihood block-support function (genotype_block_support.hpp) on synthetic groups.
#include "panvar/genotype_block_support.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool c, const std::string& m) {
    if (!c) throw std::runtime_error(m);
}

using panvar::BlockGenotypeRow;
using panvar::SupportGroup;

// 5 blocks; the called genotype is (0,1) at every block.
BlockGenotypeRow called_row() { return BlockGenotypeRow(5, {0, 1}); }
BlockGenotypeRow with(BlockGenotypeRow r, std::size_t b, std::pair<std::size_t, std::size_t> g) {
    r[b] = g;
    return r;
}
SupportGroup group(std::vector<BlockGenotypeRow> rows, double F, double E, double P, double u = 1e-9) {
    SupportGroup g;
    g.rows = std::move(rows);
    g.fragment = F;
    g.exposure = E;
    g.prior = P;
    g.bound = u;
    return g;
}

}  // namespace

int main() {
    try {
        const BlockGenotypeRow W = called_row();
        std::vector<SupportGroup> gs;
        gs.push_back(group({W}, -100.0, -10.0, -5.0));                          // 0: the call (L -110, S -115)
        gs.push_back(group({with(W, 0, {0, 2})}, -95.0, -20.0, -5.0));          // 1: wins on F, loses on L (exposure) at block 0
        gs.push_back(group({with(W, 1, {1, 1})}, -99.0, -10.0, -20.0));         // 2: wins on L at block 1; the prior selects the call
        gs.push_back(group({W}, -90.0, -10.0, -30.0));                          // 3: same genotype everywhere (phase alternative), higher L
        gs.push_back(group({with(W, 4, {2, 2})}, -100.0, -10.0 + 1e-12, -6.0)); // 4: disagrees at block 4, tied within the bounds
        const auto r = panvar::block_data_support_v1(gs, 0, 5);
        // (1) exposure: a fragment-term winner loses on the full data likelihood -> supported (legacy would say PRIOR_ONLY)
        require(r.blocks[0].label == "DATA_SUPPORTED" && r.blocks[0].certified, "gate 1: exposure-decided rival must leave block 0 DATA_SUPPORTED");
        require(std::fabs(r.blocks[0].delta_fragment + 5.0) < 1e-12 && std::fabs(r.blocks[0].delta_exposure - 10.0) < 1e-12,
                "gate 1: call-minus-competitor fragment -5, exposure +10");
        // (2) prior dependence stays visible. Group 3 (a phase alternative carrying the call's genotype with higher L) is left
        // out here: with it, block 1's genotype is data-supported through that alternative, as the frozen rule says (gate 3).
        {
            const std::vector<SupportGroup> g2{gs[0], gs[1], gs[2], gs[4]};
            const auto x = panvar::block_data_support_v1(g2, 0, 5);
            require(x.blocks[1].label == "PRIOR_DEPENDENT", "gate 2: a likelihood-winning rival selected against by the prior is PRIOR_DEPENDENT");
            require(x.blocks[1].delta_prior > 0.0 && x.blocks[1].delta_likelihood < 0.0, "gate 2: the deltas show prior for, likelihood against");
            require(r.blocks[1].label == "DATA_SUPPORTED", "gate 2/3: with the phase alternative present, block 1 is supported through it");
        }
        // (3) same genotype, competing phase: support unaffected; the phase competitor is genotype-preserving
        require(r.blocks[2].label == "NO_EVALUATED_ALTERNATIVE" && r.blocks[3].label == "NO_EVALUATED_ALTERNATIVE",
                "gate 3: a phase alternative never creates a disagreeing genotype");
        require(r.phase.genotype_preserving == 1, "gate 3: the same-genotype contender is classed genotype-preserving");
        require(r.phase.genotype_changing == 3, "gate 3: groups 1, 2 and 4 are genotype-changing fragment contenders");
        // (4) a rival disagreeing at one block affects only that block (group 1 only at block 0, group 2 only at block 1)
        require(r.blocks[2].disagreeing == 0 && r.blocks[3].disagreeing == 0, "gate 4: no disagreement leaks to other blocks");
        // (6) numerical tie inside the bounds -> DATA_TIED; the same distance outside the bounds -> certified
        // (isolated fixture)
        {
            std::vector<SupportGroup> g6{gs[0], gs[4]};   // only the call and the tied competitor: any other group carrying the
                                                          // call's genotype at block 4 with higher L would support it (frozen rule)
            const auto x = panvar::block_data_support_v1(g6, 0, 5);
            require(x.blocks[4].label == "DATA_TIED" && !x.blocks[4].certified, "gate 6: a tie within bounds is DATA_TIED");
            g6[1].exposure = -10.0 - 1e-6;   // 1e-6 below: outside 2e-9 of bounds
            require(panvar::block_data_support_v1(g6, 0, 5).blocks[4].label == "DATA_SUPPORTED", "gate 6: beyond the bounds it is certified");
        }
        // gate 3, the other direction: a same-genotype group with HIGHER L supports the genotype even when the call has lower L
        {
            std::vector<SupportGroup> g3{group({W}, -100, -10, -5), group({W}, -95, -10, -40), group({with(W, 0, {0, 3})}, -98, -10, -50)};
            const auto x = panvar::block_data_support_v1(g3, 0, 5);
            require(x.blocks[0].label == "DATA_SUPPORTED" && std::fabs(x.blocks[0].margin - 3.0) < 1e-12,
                    "gate 3: the genotype's best evaluated likelihood (a phase alternative) is the reference, margin 3");
        }
        // (5) equivalent sequences: one group with two members of one genotype row; a projection-ambiguous winner; a mixed competitor
        {
            std::vector<SupportGroup> g5{group({W, with(W, 2, {0, 2})}, -100, -10, -5),   // winner ambiguous at block 2
                                         group({W, with(W, 3, {3, 3})}, -99, -10, -30),   // mixed at block 3, not below the call
                                         group({W}, -101, -10, -6)};                      // an equivalent-genotype group
            const auto x = panvar::block_data_support_v1(g5, 0, 5);
            require(x.blocks[2].label == "AMBIGUOUS_PROJECTION", "gate 5: a projection-ambiguous winner is never resolved by a representative");
            require(x.blocks[3].label == "NO_EVALUATED_ALTERNATIVE" && x.blocks[3].mixed == 1,
                    "gate 5: a mixed group is not a disagreeing alternative");
            std::vector<SupportGroup> g6 = g5;
            g6.push_back(group({with(W, 3, {3, 4})}, -120, -10, -5));   // a definite disagreement at block 3, far below
            const auto y = panvar::block_data_support_v1(g6, 0, 5);
            require(y.blocks[3].label == "AMBIGUOUS_COMPETITOR", "gate 5: a mixed group not certifiably below blocks DATA_SUPPORTED");
        }
        // (7) order invariance: permute groups; labels, margins and the competitor's identity agree
        {
            std::vector<std::size_t> perm(gs.size());
            std::iota(perm.begin(), perm.end(), 0);
            std::reverse(perm.begin(), perm.end());
            std::vector<SupportGroup> p;
            for (const std::size_t k : perm) p.push_back(gs[k]);
            const std::size_t winner = static_cast<std::size_t>(std::find(perm.begin(), perm.end(), std::size_t{0}) - perm.begin());
            const auto q = panvar::block_data_support_v1(p, winner, 5);
            for (std::size_t b = 0; b < 5; ++b) {
                require(q.blocks[b].label == r.blocks[b].label, "gate 7: labels must not depend on group order");
                require((std::isnan(q.blocks[b].margin) && std::isnan(r.blocks[b].margin)) || q.blocks[b].margin == r.blocks[b].margin,
                        "gate 7: margins must not depend on group order");
                if (r.blocks[b].competitor >= 0) require(perm[static_cast<std::size_t>(q.blocks[b].competitor)] == static_cast<std::size_t>(r.blocks[b].competitor),
                                                         "gate 7: the strongest competitor must be the same group");
            }
            require(q.phase.genotype_changing == r.phase.genotype_changing && q.phase.genotype_preserving == r.phase.genotype_preserving,
                    "gate 7: phase classes must not depend on order");
        }
        // quality terms (frozen mapping)
        require(panvar::data_support_quality_reason("DATA_SUPPORTED").empty() && panvar::data_support_quality_reason("NO_EVALUATED_ALTERNATIVE").empty() &&
                    panvar::data_support_quality_reason("PRIOR_DEPENDENT") == "prior_dependent" && panvar::data_support_quality_reason("DATA_TIED") == "data_tied",
                "the frozen label-to-quality mapping");
        std::cout << "genotype block support: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "genotype block support: FAIL: " << e.what() << '\n';
        return 1;
    }
}
