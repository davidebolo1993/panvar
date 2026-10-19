# M3-E Amendment 1: comparator column handling and profile identity

Declared after `compare` (176.0 of 600 s spent, exit 1: 138 checks, 18 failed) and before any further execution.

## The failures

The 18 failures are 6 runs × 3 files. Neither cause is an annotation leak.

1. **`.block_genotypes.tsv`.** The comparator assumed the 13 experimental columns are the last ones and compared the first
   29 positions.
   - The headline stage appends `genotype_filter`, `genotype_filter_reasons`, `dosage_filter` and
     `dosage_filter_reasons` **after** the writer's columns, so the v1x columns sit between `quality_reason` and those
     four.
   - The positional comparison therefore compared different columns.
   - The design text said "appended", which is corrected now: the columns are added after `quality_reason`, and readers
     must use column names.
2. **`.mosaic_status.tsv` and `.settings.tsv`.** These differ only in the profile identity: the `profile` path and
   `profile_sha256`, and the settings `source` column naming the profile file. The profiles differ by design (one added
   key).

## Change

`compare.a1` makes two adjustments:
- it compares the block table **by column name**: every pre-existing column under its own name, and only `*_v1x`
  columns may be new;
- it drops `profile` and `profile_sha256` from the status comparison, and compares settings on key and value only.

Every other allowance, and the audit cross-check, is unchanged. `compare` stays recorded with exit 1.

**Unit:** `compare.a1`, declared timeout 30 s.
