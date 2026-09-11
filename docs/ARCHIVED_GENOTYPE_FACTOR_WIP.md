# Archived genotype-factor work

This checkpoint preserves the final, non-production state of the higher-order
fragment-factor experiment before development moved to the smaller mosaic
proposal-and-rescoring design.

The optimized full-domain certification fixture passes and the C4 certifier was
reduced from an estimated ~16 hours to 18.7 seconds. The audit-only serial
reporting path still has an unresolved segmentation fault at 8,000 or more
fragments, localized to the `as_frag` reporting stage. Debug stage markers are
intentionally retained in this archive.

This code is research history, not a release candidate and not the base of the
replacement genotyper. The replacement is developed from `main` on the
`genotype-mosaic` branch.
