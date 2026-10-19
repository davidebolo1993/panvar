# M1-A: freeze the performance baseline and begin failure attribution

Branch `m1-a-attribution` from `282e64e`. No production code change; no truth in any caller.

## Baseline
`baseline/` is built from the accepted records, the code and the binaries by `build_baseline.py`:
- `baseline_manifest.tsv`: code, binaries, envelope, records, explicit defaults;
- `baseline_cells.tsv`: inputs and identities;
- `baseline_settings.tsv`: every effective value with its source.

The S2-L control settings files omit the two experimental keys; their verified defaults are recorded separately.

## Audit set (`AUDIT_SET.tsv`, frozen before any new execution)
Six existing, already-evaluated cases with saved outputs:
- C4: `lzo.HG00096.5` (correct control);
- C4: `mosaic.one_homologue_once_1.2` (low-coverage mosaic error; its 5× twin is correct);
- CYP2D6: `lzo.HG02769.30` (sv:3/sv:5 labelling ambiguity, amendment-5 record);
- CYP2D6: `loo.NA19650.5` (LOO sv:9 structural error, PASS);
- ACOT: HG00097 (LZO) at 5× and 30×, with simulation truth (`truth.fa`) and no saved evaluation.

Every input exists. No campaign is rebuilt and no held-out data outside these saved evaluations is opened. The CYP2D6 LOO case is already published in the pilot record; it is development evidence from now on.

## Targets, kept separate per case
- whole diploid sequence, allowing a global homologue exchange;
- structure, where truth supports it;
- block labels and representation ambiguity;
- phase.

An availability-floor block match is never taken as whole-sequence or phase proof.

## Trace categories
`MODEL_REPRESENTATION`, `PROPOSAL`, `SEARCH`, `RANKING`, `PROJECTION`, `PHASE`, `REPORTING`, `UNDETERMINED`.

Evidence comes first from the saved records: candidates, scores, equivalence set, rounds, evaluation.

## Isolated scoring (`m1a_run.py`)
`genotype_probe score-set` (the unchanged reference scorer) runs on declared diplotypes in `candidates/`:
- every saved candidate's expected total must reproduce;
- truth is not loaded;
- for LOO, only actually scored or single-block alternatives are declared; no per-block-floor oracle is invented.

The ACOT check compares truth and called whole-sequence hashes.
