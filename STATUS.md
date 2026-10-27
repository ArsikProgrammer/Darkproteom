# Status (corrected)

## Mapping status: PARTIAL (not validated)

### Mapping table (mutations passing each filter)

| Cohort | Total SNP | CDS found | Translation match | UniProt found | Window ±5 match | In analysis |
|---|---|---|---|---|---|---|
| George | 34583 | 34408 (99.5%) | 33239 (96.6%) | 31938 (92.4%) | TBD | TBD |
| JHU | 7858 (6972 with transcript) | 0 (RefSeq) | n/a | 6900 (99.0%) | TBD | TBD |
| CLCGP | 8078 | 8066 (99.8%) | 7784 (96.5%) | 7293 (90.3%) | TBD | TBD |

### Discard fractions
- George: no UniProt 7.6%
- JHU: no UniProt 1.0%
- CLCGP: no UniProt 9.7%

### Mismatch reasons (George/CLCGP)
- To be classified: complementary strand, codon phase shift, transcript mismatch.

### Duplicate samples
- George vs CLCGP: S01512 (overlap 0.989) → keep George, remove CLCGP.
- JHU vs CLCGP: H1672 (overlap 0.675) → manual check required.

### SBS4 stability
- CI width = 0.19 < 0.3 → not unstable.
- Whole-genome vs coding: NOT performed (George "all" ≈ coding, comparison invalid).

### GEE model
- Synthetic: OR=1.5 → 1.75 (single run). Need 500+ runs for bias/coverage.
- Status: OPEN.

### Power simulation
- Still overestimates (sample-level, not mutation-level).
- Status: OPEN.

### Archive
- Download in progress (FTP, resume).
- Status: IN PROGRESS.
