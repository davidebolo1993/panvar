#include "panvar/genotype_joint_search.hpp"

#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar {
namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point mark) { return std::chrono::duration<double>(Clock::now() - mark).count(); }

std::size_t shortlisted(const PreparedGenotypeModel& model, std::size_t block, std::size_t top) {
    return std::min(top, model.blocks()[block].alleles.size() - 1);
}

void add_work(JointSearchResult& r, const BlockDeltaWork& w) {
    r.origins_reused += w.origins_reused;
    r.origins_recomputed += w.origins_recomputed;
    r.dp_executions += w.dp_executions;
    r.window_bases += w.local_window_bases;
}

// The gate hook's deterministic permutation (a fixed LCG; production never sets a seed).
template <typename T>
void permute(std::vector<T>& items, std::uint64_t seed) {
    if (seed == 0) return;
    std::uint64_t s = seed;
    for (std::size_t i = items.size(); i > 1; --i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        std::swap(items[i - 1], items[(s >> 33) % i]);
    }
}

// Issues the requests in the given order and returns the slots in the ORIGINAL order.
std::vector<std::size_t> issue(PreparedScoringContext& context, std::vector<RepairRequest> requests, std::uint64_t seed,
                               JointSearchResult& result) {
    std::vector<std::size_t> order(requests.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    permute(order, seed);
    std::vector<RepairRequest> issued;
    issued.reserve(requests.size());
    for (const std::size_t i : order) issued.push_back(requests[i]);
    std::vector<BlockDeltaWork> works;
    const std::vector<std::size_t> got = context.add_repaired_haplotypes(issued, &works);
    for (const BlockDeltaWork& w : works) add_work(result, w);
    std::vector<std::size_t> out(requests.size());
    for (std::size_t k = 0; k < order.size(); ++k) out[order[k]] = got[k];
    return out;
}

}  // namespace

JointSearchBounds joint_search_bounds(const PreparedGenotypeModel& model, const JointSearchOptions& options) {
    JointSearchBounds b;
    const std::size_t blocks = model.blocks().size();
    for (std::size_t i = 0; i < blocks; ++i) b.single_alternatives += 2 * (model.blocks()[i].alleles.size() - 1);
    for (int h = 0; h < 2; ++h) {
        for (std::size_t i = 0; i < blocks; ++i) b.catalogue_raw += shortlisted(model, i, options.top);
        for (std::size_t i = 0; i < blocks; ++i) {
            for (std::size_t j = i + 1; j < blocks; ++j) {
                b.catalogue_raw += static_cast<std::uint64_t>(shortlisted(model, i, options.top)) * shortlisted(model, j, options.top);
            }
        }
    }
    b.within = b.single_alternatives <= options.single_limit && b.catalogue_raw <= options.catalogue_limit;
    return b;
}

JointSearchResult run_joint_block_search(
    PreparedScoringContext& context, const PreparedGenotypeModel& model, const std::vector<std::size_t>& first,
    const std::vector<std::size_t>& second, std::size_t slot_first, std::size_t slot_second, double incumbent_total,
    std::size_t fragments, std::size_t mismatch_band, const JointSearchOptions& options,
    const std::function<void(const JointCandidate&, const std::string&, const std::string&)>& on_scored,
    const std::vector<std::size_t>& keep) {
    const Clock::time_point start = Clock::now();
    JointSearchResult r;
    r.incumbent_total = incumbent_total;
    r.bounds = joint_search_bounds(model, options);
    if (!r.bounds.within) {
        r.status = "skipped_bounds";
        r.decision = "skipped";
        return r;
    }
    const std::size_t blocks = model.blocks().size();
    if (first.size() != blocks || second.size() != blocks) throw std::invalid_argument("joint search: incumbent paths and model differ");
    if (!context.holds_census(slot_first) || !context.holds_census(slot_second)) {
        throw std::logic_error("joint search: the incumbent slots hold no census");
    }
    const std::vector<std::size_t> inc[2] = {first, second};
    const std::size_t islot[2] = {slot_first, slot_second};
    const HaplotypeSequence spelled[2] = {spell_mosaic(model, MosaicPath{first}), spell_mosaic(model, MosaicPath{second})};
    const std::uint64_t placements_before = context.full_placement_count();
    const std::vector<double> switches(blocks - 1, options.switch_probability);

    // ---- 1. every single-block alternative, other homologue fixed (shortlisting only)
    Clock::time_point mark = Clock::now();
    struct Single {
        int h;
        std::size_t b, a, slot = 0, rank = 0;
        double likelihood = 0.0;
    };
    std::vector<Single> singles;
    std::vector<RepairRequest> requests;
    for (int h = 0; h < 2; ++h) {
        for (std::size_t b = 0; b < blocks; ++b) {
            const SequenceInterval iv = spelled[h].block_intervals[b];
            for (std::size_t a = 0; a < model.blocks()[b].alleles.size(); ++a) {
                if (a == inc[h][b]) continue;
                singles.push_back({h, b, a});
                requests.push_back({islot[h], iv.begin, iv.end, model.blocks()[b].alleles[a].sequence});
            }
        }
    }
    {
        const std::vector<std::size_t> got = issue(context, requests, options.request_order_seed, r);
        std::vector<HaplotypeCandidate> pairs;
        for (std::size_t i = 0; i < singles.size(); ++i) {
            singles[i].slot = got[i];
            pairs.push_back(singles[i].h == 0 ? HaplotypeCandidate{got[i], islot[1]} : HaplotypeCandidate{islot[0], got[i]});
        }
        const FragmentScoringResult scored = context.score_candidates(pairs);
        for (std::size_t i = 0; i < singles.size(); ++i) singles[i].likelihood = scored.scores[i].log_likelihood;
    }
    r.singles_scored = singles.size();
    r.seconds_singles = since(mark);
    std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> shortlist;
    {
        std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> group;
        for (std::size_t i = 0; i < singles.size(); ++i) group[{singles[i].h, singles[i].b}].push_back(i);
        for (auto& [key, idx] : group) {
            std::sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) {
                return singles[x].likelihood != singles[y].likelihood ? singles[x].likelihood > singles[y].likelihood
                                                                      : singles[x].a < singles[y].a;
            });
            for (std::size_t k = 0; k < idx.size(); ++k) {
                singles[idx[k]].rank = k + 1;
                if (k < options.top) shortlist[key].push_back(idx[k]);
            }
        }
    }

    // ---- 2. promote the shortlisted singles: census by delta (materialised), as the repair promotes a winner
    mark = Clock::now();
    for (const auto& [key, idx] : shortlist) {
        const SequenceInterval iv = spelled[key.first].block_intervals[key.second];
        for (const std::size_t i : idx) {
            BlockDeltaWork w;
            const std::size_t slot = context.add_repaired_haplotype(islot[key.first], iv.begin, iv.end,
                                                                    model.blocks()[key.second].alleles[singles[i].a].sequence, &w, true);
            if (slot != singles[i].slot) throw std::logic_error("joint search: promotion did not return the scored slot");
            add_work(r, w);
            ++r.promoted;
        }
    }
    r.census_bytes_peak = context.retained_census_bytes();
    r.seconds_promotion = since(mark);

    // ---- 3. the catalogue: shortlisted singles, then every same-homologue block pair from the shortlists
    mark = Clock::now();
    const auto move_text = [&](int h, std::size_t b, std::size_t from, std::size_t to, std::size_t rank) {
        std::ostringstream m;
        m << "h" << h << ":" << model.blocks()[b].block_id << ":" << from << ">" << to << ":r" << rank;
        return m.str();
    };
    struct Pending {
        int h;
        std::vector<std::size_t> path;
        std::string move;
        std::size_t slot;
    };
    std::vector<Pending> pending;
    for (const auto& [key, idx] : shortlist) {
        for (const std::size_t i : idx) {
            std::vector<std::size_t> p = inc[key.first];
            p[key.second] = singles[i].a;
            pending.push_back({key.first, p, "single:" + move_text(key.first, key.second, inc[key.first][key.second], singles[i].a, singles[i].rank),
                               singles[i].slot});
        }
    }
    std::vector<RepairRequest> joint_requests;
    std::vector<Pending> joint_pending;
    for (int h = 0; h < 2; ++h) {
        for (std::size_t b = 0; b < blocks; ++b) {
            for (std::size_t c = b + 1; c < blocks; ++c) {
                const auto kb = shortlist.find({h, b}), kc = shortlist.find({h, c});
                if (kb == shortlist.end() || kc == shortlist.end()) continue;
                const SequenceInterval ic = spelled[h].block_intervals[c];
                for (const std::size_t i : kb->second) {
                    const long shift = static_cast<long>(model.blocks()[b].alleles[singles[i].a].sequence.size()) -
                                       static_cast<long>(model.blocks()[b].alleles[inc[h][b]].sequence.size());
                    for (const std::size_t j : kc->second) {
                        joint_requests.push_back({singles[i].slot, static_cast<std::size_t>(static_cast<long>(ic.begin) + shift),
                                                  static_cast<std::size_t>(static_cast<long>(ic.end) + shift),
                                                  model.blocks()[c].alleles[singles[j].a].sequence});
                        std::vector<std::size_t> p = inc[h];
                        p[b] = singles[i].a;
                        p[c] = singles[j].a;
                        joint_pending.push_back({h, p, "pair:" + move_text(h, b, inc[h][b], singles[i].a, singles[i].rank) + "+" +
                                                           move_text(h, c, inc[h][c], singles[j].a, singles[j].rank).substr(3), 0});
                    }
                }
            }
        }
    }
    r.joint_requests = joint_requests.size();
    {
        const std::vector<std::size_t> got = issue(context, joint_requests, options.request_order_seed, r);
        for (std::size_t k = 0; k < joint_pending.size(); ++k) {
            joint_pending[k].slot = got[k];
            pending.push_back(joint_pending[k]);
        }
    }
    r.transient_bytes_peak = context.transient_payload_bytes();
    r.seconds_joint = since(mark);
    // deduplicate by unordered slot pair: slots are byte-confirmed sequence identities
    const auto pair_of = [&](const Pending& p) {
        const std::size_t a = p.h == 0 ? p.slot : islot[0], b = p.h == 1 ? p.slot : islot[1];
        return std::make_pair(a, b);
    };
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> seen;
    const auto incumbent_key = std::minmax(islot[0], islot[1]);
    for (const Pending& p : pending) {
        const auto slots = pair_of(p);
        const auto key = std::minmax(slots.first, slots.second);
        if (key == incumbent_key) continue;
        const auto found = seen.find(key);
        if (found != seen.end()) {
            ++r.catalogue[found->second].equivalent_moves;
            continue;
        }
        seen.emplace(key, r.catalogue.size());
        JointCandidate c;
        c.first = p.h == 0 ? p.path : inc[0];
        c.second = p.h == 1 ? p.path : inc[1];
        c.slot_first = slots.first;
        c.slot_second = slots.second;
        c.move = p.move;
        r.catalogue.push_back(std::move(c));
    }

    // ---- 4. complete scores: combination from retained columns + the sequence-level pair prior
    mark = Clock::now();
    std::vector<HaplotypeCandidate> pairs{{std::min(islot[0], islot[1]), std::max(islot[0], islot[1])}};
    for (const JointCandidate& c : r.catalogue) pairs.push_back({c.slot_first, c.slot_second});
    const FragmentScoringResult scored = context.score_candidates(pairs);
    r.seconds_scoring = since(mark);
    mark = Clock::now();
    std::map<std::vector<std::size_t>, std::pair<double, std::string>> haploid;   // path -> (sequence prior, bytes)
    const auto haploid_of = [&](const std::vector<std::size_t>& p) -> const std::pair<double, std::string>& {
        auto f = haploid.find(p);
        if (f == haploid.end()) {
            std::string bytes = spell_mosaic(model, MosaicPath{p}).bases;
            const double prior = log_haploid_sequence_prior(model, bytes, switches);
            f = haploid.emplace(p, std::make_pair(prior, std::move(bytes))).first;
        }
        return f->second;
    };
    const auto pair_prior = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
        const auto& x = haploid_of(a);
        const auto& y = haploid_of(b);
        if (!std::isfinite(x.first) || !std::isfinite(y.first)) return -std::numeric_limits<double>::infinity();
        return x.first + y.first + (x.second != y.second ? std::log(2.0) : 0.0);   // log 2 decided on bytes
    };
    std::vector<DiplotypeScore> scores{scored.scores[0]};
    std::vector<double> priors{pair_prior(first, second)};
    for (std::size_t i = 0; i < r.catalogue.size(); ++i) {
        r.catalogue[i].score = scored.scores[i + 1];
        r.catalogue[i].prior = pair_prior(r.catalogue[i].first, r.catalogue[i].second);
        scores.push_back(r.catalogue[i].score);
        priors.push_back(r.catalogue[i].prior);
    }
    r.seconds_priors = since(mark);
    // the incumbent must reproduce (the trunk's declared total) before any move is interpreted
    const double reproduced = scores[0].log_likelihood + priors[0];
    if (reproduced != incumbent_total) {   // explicit refusal: nothing is moved or reported, every pass slot is released
        r.status = "refused";
        r.decision = "incumbent_not_reproduced";
        std::vector<std::size_t> kept = keep;
        kept.push_back(islot[0]);
        kept.push_back(islot[1]);
        r.slots_created = context.live_transient_slots().size();
        context.release_transient(kept);
        r.transient_bytes_after = context.transient_payload_bytes();
        r.live_transient_after = context.live_transient_slots().size();
        r.full_placements = context.full_placement_count() - placements_before;
        r.seconds_total = since(start);
        return r;
    }

    // ---- 5. the repair round's certified move rule, over the incumbent and the catalogue
    const CertifiedRanking ranking =
        certify_scored_ranking(scores, fragments, mismatch_band, scored.counters, CertificationOptions{}, priors);
    const std::size_t top = ranking.candidates.front().candidate_id;
    r.winner_set_size = ranking.winner_set.size();
    r.certified_gap = ranking.certified_gap;
    for (const std::size_t id : ranking.winner_set) r.incumbent_in_top_set = r.incumbent_in_top_set || id == 0;
    if (ranking.status != RescoreStatus::Complete) {
        r.decision = "rescore_incomplete";
    } else if (top == 0) {
        r.decision = "no_improvement";
    } else if (ranking.winner_set.size() != 1) {
        r.decision = "tie_at_best";
    } else if (scores[top].log_likelihood + priors[top] <= incumbent_total) {
        r.decision = "no_improvement";
    } else {
        r.decision = "moved";
        r.moved = true;
        r.winner = top - 1;
    }
    r.status = "complete";

    // ---- 6. report every catalogue member (complete ranking only), then release
    if (on_scored) {
        for (const JointCandidate& c : r.catalogue) on_scored(c, haploid_of(c.first).second, haploid_of(c.second).second);
    }
    std::vector<std::size_t> kept = keep;
    kept.push_back(islot[0]);
    kept.push_back(islot[1]);
    if (r.moved) {
        // The winner keeps a census for the continuation. A joint winner was scored without one: retain it now, by the
        // same delta from its promoted single (fills the existing slot; no full placement).
        const JointCandidate& w = r.catalogue[r.winner];
        for (const std::size_t s : {w.slot_first, w.slot_second}) {
            if (context.holds_census(s)) continue;
            bool done = false;
            for (std::size_t k = 0; k < joint_pending.size() && !done; ++k) {
                if (joint_pending[k].slot != s) continue;
                (void)context.add_repaired_haplotype(joint_requests[k].incumbent, joint_requests[k].block_begin,
                                                     joint_requests[k].block_end, joint_requests[k].replacement, nullptr, true);
                done = true;
            }
            if (!done || !context.holds_census(s)) throw std::logic_error("joint search: the winner's census could not be retained");
        }
        kept.push_back(w.slot_first);
        kept.push_back(w.slot_second);
    }
    r.slots_created = context.live_transient_slots().size();
    context.release_transient(kept);
    r.transient_bytes_after = context.transient_payload_bytes();
    r.live_transient_after = context.live_transient_slots().size();
    r.full_placements = context.full_placement_count() - placements_before;
    r.seconds_total = since(start);
    return r;
}

}  // namespace panvar
