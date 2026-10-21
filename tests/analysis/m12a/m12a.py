#!/usr/bin/env python3
"""M12-A: repeated-unit placement/exposure consistency audit (diagnostic; truth used only for diagnostic checks).

Modes:
  contract <out>   rate/exposure provenance per case: declared library vs simulator/subsample provenance; marker rate and
                   count depth against the derived expectation; locus-level observed vs expected fragments
  prepare <out>    per case: score-set candidates (incumbent with expected total, assign_A, explicit truth-in-context in both
                   homologue assignments) and count-attrib candidates (incumbent, best_A, truth block in both assignments)
  analyse <out>    likelihood and count comparisons, per-fragment scale split, block-level observed vs expected fragments,
                   and the frozen classification (DESIGN.md §5)
Reads the sealed M11-D/M11-E/M6-A packages; writes nothing there.
"""
import functools
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m11e  # noqa: E402

m11d = m11e.m11d
T, KV, write = m11d.T, m11d.KV, m11d.write
fasta = functools.lru_cache(maxsize=None)(m11d.fasta)
R_ = m11d.R
M11D, M11E = m11e.M11D, m11e.OUT
OUT = f"{R_}/m12a_rate_exposure_audit_20261021"
K = 31                     # marker / count k-mer length (caller marker panel k=31,s=15; M6-A count k=31)
LN2 = math.log(2.0)        # per-fragment scale split: a copy-count change moves a fragment by at most ln 2 (k >= 1)
RATE_TOL = 0.05            # rate-contract tolerance (relative)
Z_LOCUS = 4.0              # locus-level count tolerance (Poisson z)
Z_BLOCK = 3.0              # block-level count: a sampling-compatible deviation

# frozen cases (DESIGN.md §2): M11-E K1, K4, K5 and its three controls; K4's count evidence is M6-A's (reused, not rerun)
CASES = [dict(case="K1", cell="c4.NA24385", block="sv:7", role="catastrophic", count="run"),
         dict(case="K4", cell="cyp.NA19650", block="sv:9", role="catastrophic", count="m6a"),
         dict(case="K5", cell="cyp.NA19983", block="sv:9", role="catastrophic", count="run"),
         dict(case="R1", cell="cyp.HG02953", block="sv:9", role="control", count="run"),
         dict(case="R2", cell="c4.HG00358", block="sv:7", role="control", count="run"),
         dict(case="R3", cell="acot.HG03209", block="sv:13", role="control", count="run")]
M6A_K4 = f"{R_}/m6a_count_model_20261019/runs/countfit.cyp_loo_na19650_5/cf"


def case(cid):
    return next(c for c in CASES if c["case"] == cid)


def provenance(cell):
    """Simulation directory (simulation.tsv, truth.fa) and the subsample record of the reads the run used."""
    c = m11d.cell(cell)
    d = c["donor"]
    reads = m11d.inputs(c)["reads"]
    if c["locus"] == "C4":
        sim = f"{R_}/genotype_coverage_c4_held_out_abb35a0/{d}/source"
    elif c["locus"] == "CYP2D6":
        sim = f"{R_}/cyp2d6_shortread_pilot_v1/data/reads/{d}/cov30"
    else:
        sim = f"{R_}/m8c_copying_challenge_20261020/runs/prepare.a1/reads/{d}/cov30"
    return dict(reads=reads, sim=sim, subsample=f"{reads}/subsample.tsv", library=f"{reads}/library.txt")


def library(p):
    return {l.split()[0]: float(l.split()[1]) for l in open(p) if l.strip() and not l.startswith("#")}


def read_stats(reads):
    n, bases, mates = 0, 0, 0
    for name in ("reads_1.fq", "reads_2.fq"):
        with open(f"{reads}/{name}") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    bases += len(line.rstrip("\n"))
                    mates += 1
        n = mates // 2
    return n, bases / mates


def insert_mean(lib):
    """Mean of the declared discretised insert distribution (normal, mean +- sigmas*sd, integer support)."""
    mu, sd, k = lib["insert_mean_bp"], lib["insert_sd_bp"], lib["support_sigmas"]
    lo, hi = int(math.ceil(mu - k * sd)), int(math.floor(mu + k * sd))
    w = [(i, math.exp(-0.5 * ((i - mu) / sd) ** 2)) for i in range(lo, hi + 1)]
    return sum(i * x for i, x in w) / sum(x for _, x in w)


def expected_kmer_rate(lib, read_len, k=K):
    """Expected exact observations of one k-mer per chromosome copy, both mates, canonical (strand-independent):
    2 mates x (1-eta) lambda templates/base/copy x (R-k+1) covering starts x (1-e)^k error-free (interior; boundaries apart)."""
    return 2.0 * (1.0 - lib["outlier_weight"]) * lib["fragment_start_intensity"] * (read_len - k + 1) * (1.0 - lib["substitution_error"]) ** k


# ------------------------------------------------------------------------------------------------ contract
def contract(out):
    rows = []
    for c in CASES:
        p = provenance(c["cell"])
        lib, sim, sub = library(p["library"]), KV(f"{p['sim']}/simulation.tsv"), KV(p["subsample"])
        n, rlen = read_stats(p["reads"])
        st = KV(f"{M11D}/runs/call.{c['cell']}/call.mosaic_status.tsv")
        frac = float(sub["requested_fraction"])
        lam_sim = float(sim["lambda"])
        exp_templates_30 = float(sim["haplotype_0_expected_templates"]) + float(sim["haplotype_1_expected_templates"])
        exp5 = exp_templates_30 * frac
        src = float(sub.get("source_fragments", sim.get("fragments_written", "nan")))
        kept_exp = src * frac
        rate_exp = expected_kmer_rate(lib, rlen)
        rate_obs = float(st["marker_depth_rate_per_copy"])
        rows.append(dict(
            case=c["case"], cell=c["cell"], lambda_declared=lib["fragment_start_intensity"], lambda_simulated=lam_sim,
            requested_fraction=frac, lambda_declared_over_sim_x_fraction=lib["fragment_start_intensity"] / (lam_sim * frac),
            read_length_reads=rlen, read_length_simulated=sim.get("read_length", "NA"), depth_simulated=sim.get("depth", "NA"),
            substitution_error=lib["substitution_error"], outlier_weight=lib["outlier_weight"],
            insert_mean_discrete=round(insert_mean(lib), 6),
            simulator_signal_rate_per_copy=lam_sim * (1 - lib["outlier_weight"]),
            simulator_check_depth_over_4R=float(sim.get("depth", "nan")) / (4.0 * rlen),
            truth_bp=f"{sim['haplotype_0_bp']}+{sim['haplotype_1_bp']}",
            simulated_templates_30x=f"{sim['haplotype_0_templates']}+{sim['haplotype_1_templates']}",
            expected_templates_30x=round(exp_templates_30, 1),
            fragments_observed=n, fragments_expected_truth=round(exp5, 1), z_vs_truth=round((n - exp5) / math.sqrt(exp5), 2),
            subsample_kept_expected=round(kept_exp, 1), z_subsample=round((n - kept_exp) / math.sqrt(kept_exp * (1 - frac)), 2),
            marker_rate_status=rate_obs, marker_rate_expected=round(rate_exp, 6), marker_rate_ratio=round(rate_obs / rate_exp, 4),
            marker_rate_ratio_at_kept_fraction=round(rate_obs / (rate_exp * n / kept_exp), 4)))
    write(out, "contract", rows)
    for r in rows:
        print(r["case"], r["cell"], "lambda decl/(sim*frac)", round(r["lambda_declared_over_sim_x_fraction"], 9), "rate ratio",
              r["marker_rate_ratio"], "z truth", r["z_vs_truth"])
    return 0


# ------------------------------------------------------------------------------------------------ prepare
def prepare(out):
    for c in CASES:
        cid, cell, bid = c["case"], c["cell"], c["block"]
        ev = m11e.eval_rows(cell)
        ids = [e["block_id"] for e in ev]
        b = ids.index(bid)
        win = next(r for r in T(f"{M11D}/runs/call.{cell}/call.equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true")
        h1, h2 = m11e.vec(win["hap1_alleles"]), m11e.vec(win["hap2_alleles"])
        e = ev[b]
        b1, b2 = int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"])
        seg1, seg2 = m11e.spell(cell, ids, h1), m11e.spell(cell, ids, h2)
        tr = fasta(f"{M11D}/runs/truth.{cell}/bs.{bid}.truth.fa")
        t1, t2 = tr["truth1"], tr["truth2"]
        # explicit truth-in-context: only the block's bytes replaced by the truth's (two assignments)
        recs, lines = {}, []
        for lab, x, y in (("truth_A", t1, t2), ("truth_B", t2, t1)):
            p1 = "".join(seg1[:b]) + x + "".join(seg1[b + 1:])
            p2 = "".join(seg2[:b]) + y + "".join(seg2[b + 1:])
            recs[f"{lab}.1"], recs[f"{lab}.2"] = p1, p2
            o1, o2 = len("".join(seg1[:b])), len("".join(seg2[:b]))
            lines.append(f"{lab}\t{out}/{cid}.explicit.fa\t{lab}.1\t{lab}.2\t{o1}\t{o1 + len(x)}\t{o2}\t{o2 + len(y)}\n")
        with open(f"{out}/{cid}.explicit.fa", "w") as f:
            for k, s in recs.items():
                f.write(f">{k}\n{s}\n")
        open(f"{out}/{cid}.explicit.tsv", "w").writelines(lines)
        a1, a2 = list(h1), list(h2)
        a1[b], a2[b] = b1, b2
        with open(f"{out}/{cid}.candidates.tsv", "w") as f:
            f.write(f"incumbent\t{','.join(map(str, h1))}\t{','.join(map(str, h2))}\t{win['score']}\n")
            f.write(f"assign_A\t{','.join(map(str, a1))}\t{','.join(map(str, a2))}\tNA\n")
        # count-attrib: alleles in the run's context; the truth block bytes (oracle) in both assignments
        with open(f"{out}/{cid}.block.fa", "w") as f:
            f.write(f">truth1\n{t1}\n>truth2\n{t2}\n")
        with open(f"{out}/{cid}.count_candidates.tsv", "w") as f:
            f.write(f"incumbent\talleles\t{h1[b]}\t{h2[b]}\t0\nbest_A\talleles\t{b1}\t{b2}\t0\n"
                    f"truth_A\tblock\ttruth1\ttruth2\t1\ntruth_B\tblock\ttruth2\ttruth1\t1\n")
        sim = KV(f"{provenance(cell)['sim']}/simulation.tsv")
        tf = fasta(f"{provenance(cell)['sim']}/truth.fa")
        with open(f"{out}/{cid}.truth_full.fa", "w") as f:
            f.write(f">truth.1\n{tf[sim['haplotype_0']]}\n>truth.2\n{tf[sim['haplotype_1']]}\n")
        print(cid, cell, bid, "block alleles", (h1[b], h2[b]), "best", (b1, b2), "truth block bp", (len(t1), len(t2)))
    return 0


# ------------------------------------------------------------------------------------------------ analyse
def scale_split(path):
    small = large = 0.0
    n_small = n_large = 0
    for r in T(path):
        d = float(r["delta"])
        if abs(d) <= LN2 + 1e-9:
            small += d
            n_small += d != 0.0
        else:
            large += d
            n_large += 1
    return small, n_small, large, n_large


def block_counts(summary, side):
    """Fragments with any placement touching the block on candidate side A or B (decomposition classes)."""
    g = lambda k: float(summary.get(k, 0) or 0)  # noqa: E731
    only = g(f"class_block_on_{side}_only_fragments")
    return only + g("class_block_interior_on_both_fragments") + g("class_block_junction_on_both_fragments")


def analyse(out):
    contract_rows = {r["case"]: r for r in T(f"{OUT}/runs/contract/contract.tsv")}
    rows = []
    fixtures_ok = os.path.exists(f"{OUT}/runs/gate.fixtures/run.log") and "ALL GATES PASS" in open(f"{OUT}/runs/gate.fixtures/run.log").read()
    for c in CASES:
        cid, cell = c["case"], c["cell"]
        ss = f"{OUT}/runs/lik.{cid}/ss"
        if not os.path.exists(ss + ".scores_truth_free.tsv"):
            rows.append(dict(case=cid, cell=cell, classification="UNRESOLVED_CAUSE", reason="no likelihood output"))
            continue
        sc = {r["label"]: r for r in T(ss + ".scores_truth_free.tsv")}
        lib = library(provenance(cell)["library"])
        rate = (1 - lib["outlier_weight"]) * lib["fragment_start_intensity"]
        mu = insert_mean(lib)
        con = contract_rows[cid]
        n_obs = float(con["fragments_observed"])
        inc = sc["incumbent"]
        lik = lambda r: float(r["fragment_term"]) + float(r["exposure_term"])  # noqa: E731
        reproduced = abs(float(inc["total"]) - float(inc["expected_total"])) <= 1e-6
        row = dict(case=cid, cell=cell, block_id=c["block"], role=c["role"], reproduced=reproduced, fixtures_pass=fixtures_ok,
                   rate_contract_ratio=con["marker_rate_ratio"], lambda_contract=con["lambda_declared_over_sim_x_fraction"],
                   locus_z_vs_truth=con["z_vs_truth"])
        for lab in ("incumbent", "assign_A", "truth_A", "truth_B"):
            r = sc[lab]
            row[f"{lab}_expected_fragments_locus"] = round(-float(r["exposure_term"]), 2)
            row[f"{lab}_locus_z"] = round((n_obs + 0.0 - (-float(r["exposure_term"]))) / math.sqrt(-float(r["exposure_term"])), 2)
            if lab != "incumbent":
                row[f"{lab}_d_fragment"] = float(r["fragment_term"]) - float(inc["fragment_term"])
                row[f"{lab}_d_exposure"] = float(r["exposure_term"]) - float(inc["exposure_term"])
                row[f"{lab}_d_likelihood"] = lik(r) - lik(inc)
        best_truth = max(("truth_A", "truth_B"), key=lambda l: lik(sc[l]))
        row["truth_ctx_better_assignment"] = best_truth
        row["truth_ctx_d_likelihood"] = lik(sc[best_truth]) - lik(inc)
        recon = True
        for lab in ("assign_A", best_truth):
            dp = f"{ss}.decomp.{lab}__incumbent"
            s = KV(dp + ".summary.tsv")
            ok = s.get("placement_evidence_interpretable") == "true" and all(
                abs(float(s[k])) <= 1e-6 for k in ("fragment_term_A_from_canonical_minus_scorer", "fragment_term_B_from_canonical_minus_scorer"))
            recon = recon and ok
            sm, ns, lg, nl = scale_split(dp + ".fragments.tsv")
            tag = "assign_A" if lab == "assign_A" else "truth_ctx"
            row[f"{tag}_d_fragment_multiplicity_scale"] = sm
            row[f"{tag}_n_multiplicity_scale"] = ns
            row[f"{tag}_d_fragment_emission_scale"] = lg
            row[f"{tag}_n_emission_scale"] = nl
            bp_alt = float(s["A1_block_bp"]) + float(s["A2_block_bp"])
            bp_inc = float(s["B1_block_bp"]) + float(s["B2_block_bp"])
            obs_alt, obs_inc = block_counts(s, "A"), block_counts(s, "B")
            e_alt, e_inc = rate * (bp_alt + 2 * (mu - 1)), rate * (bp_inc + 2 * (mu - 1))
            row[f"{tag}_block_bp"] = bp_alt
            row["incumbent_block_bp"] = bp_inc
            row[f"{tag}_block_fragments_observed"] = obs_alt
            row[f"{tag}_block_fragments_expected"] = round(e_alt, 1)
            row[f"{tag}_block_z"] = round((obs_alt - e_alt) / math.sqrt(e_alt), 2)
            row[f"incumbent_block_fragments_observed_vs_{tag}"] = obs_inc
            row["incumbent_block_fragments_expected"] = round(e_inc, 1)
            row[f"incumbent_block_z_vs_{tag}"] = round((obs_inc - e_inc) / math.sqrt(e_inc), 2)
        row["decompositions_reconciled"] = recon
        # count evidence (reported beside the likelihood, never added to it)
        if c["count"] == "m6a":
            cs = KV(M6A_K4 + ".summary.tsv")
            cand = T(M6A_K4 + ".candidates.tsv")
            best_a = next((r for r in cand if (r["allele_h1"], r["allele_h2"]) == tuple(row_pair(cid, "best"))), None)
            caller = next((r for r in cand if r["is_caller"] == "true"), None)
            row.update(count_source="M6-A count-fit (reused)", count_depth=cs["depth_haploid_kmer_mean"], count_stable_kmers=cs["stable_copy_kmers"],
                       count_d_loss_best_A=float(best_a["loss"]) - float(caller["loss"]) if best_a and caller else "NA",
                       count_d_loss_truth_ctx="NOT_MEASURED_IN_M6A", count_truth_multiplicity="NOT_MEASURED_IN_M6A")
        else:
            cp = f"{OUT}/runs/count.{cid}/ca"
            if os.path.exists(cp + ".candidates.tsv"):
                cs = KV(cp + ".summary.tsv")
                cand = {r["label"]: r for r in T(cp + ".candidates.tsv")}
                profiles = [k[5:] for k in next(iter(cand.values())) if k.startswith("loss_")]
                row.update(count_source="count-attrib (this batch)", count_depth=cs.get("depth_haploid_kmer_mean"),
                           count_stable_kmers=cs.get("stable_copy_kmers"),
                           count_truth_multiplicity=cs.get("stable_kmers_truth_mean_multiplicity", "NA"),
                           count_truth_kmers_outside_dictionary=cs.get("truth_kmers_outside_dictionary", "NA"))
                for pr in profiles:
                    li = float(cand["incumbent"][f"loss_{pr}"])
                    row[f"count_{pr}_d_loss_best_A"] = float(cand["best_A"][f"loss_{pr}"]) - li
                    row[f"count_{pr}_d_loss_truth_ctx"] = min(float(cand["truth_A"][f"loss_{pr}"]), float(cand["truth_B"][f"loss_{pr}"])) - li
                for x in T(cp + ".contrasts.tsv"):
                    if x["profile"] == profiles[0]:
                        key = "best_A" if x["contrast"].startswith("best_A") else "truth"
                        row[f"count_contrast_{key}_kmers_differing"] = x["kmers_differing"]
                        row[f"count_contrast_{key}_{x['category']}"] = f"{x['category_kmers']}:{x['category_delta']}"
            else:
                row["count_source"] = "NOT_MEASURED (count unit did not complete)"
        if row.get("count_depth") not in (None, "NA"):
            m = row.get("count_truth_multiplicity")
            corr = float(m) / 2.0 if m not in (None, "NA", "NOT_MEASURED_IN_M6A") else 1.0
            row["count_depth_expected"] = round(expected_kmer_rate(lib, float(con["read_length_reads"])) * corr, 6)
            row["count_depth_ratio"] = round(float(row["count_depth"]) / row["count_depth_expected"], 4)
        row["classification"], row["reason"] = classify(row)
        rows.append(row)
        print(cid, cell, row["classification"], f"truth_ctx dLik {row['truth_ctx_d_likelihood']:+.2f}", f"assign_A dLik {row['assign_A_d_likelihood']:+.2f}")
    write(out, "audit", rows)
    return 0


def row_pair(cid, which):
    c = case(cid)
    e = next(x for x in m11e.eval_rows(c["cell"]) if x["block_id"] == c["block"])
    return [e["truth1_closest_allele"], e["truth2_closest_allele"]]


def classify(r):
    """Frozen order (DESIGN.md §5)."""
    if not (r["fixtures_pass"] and r["reproduced"] and r["decompositions_reconciled"]):
        return "IMPLEMENTATION_OR_NORMALIZATION_INCONSISTENCY", "a fixture gate, the reproduction, or an origin re-summation failed"
    rate_bad = abs(float(r["rate_contract_ratio"]) - 1) > RATE_TOL or abs(float(r["lambda_contract"]) - 1) > 1e-6 \
        or abs(float(r["locus_z_vs_truth"])) > Z_LOCUS or ("count_depth_ratio" in r and abs(r["count_depth_ratio"] - 1) > RATE_TOL)
    if rate_bad:
        return "INCONSISTENT_RATE_OR_CONVENTION", "declared/simulated intensity, marker rate, count depth or locus count outside tolerance"
    if r["role"] == "control":
        return "CONTROL", "rates and reproduction checked; truth-in-context contrast reported"
    if r["truth_ctx_d_likelihood"] > 0:
        return "SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH", \
            "with the truth's block content in the winner's context the likelihood prefers the true copy number; the available alleles' content decides"
    z = r["truth_ctx_block_z"]
    toward_incumbent = (r["incumbent_block_bp"] > r["truth_ctx_block_bp"]) == (z > 0)
    if abs(z) <= Z_BLOCK and toward_incumbent:
        return "SAMPLING_COMPATIBLE_DISAGREEMENT", "even the truth content loses; the block's fragment count deviates toward the incumbent within sampling"
    if not toward_incumbent:
        return "SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH", \
            "the block's fragment count favours the true length, yet the placement likelihood prefers the incumbent with true content"
    return "UNRESOLVED_CAUSE", "truth content loses and the block count deviation exceeds the sampling tolerance"


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"contract": contract, "prepare": prepare, "analyse": analyse}[mode](out))
