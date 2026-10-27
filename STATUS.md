# Status (corrected)

## Mapping status: PARTIAL (not validated)

### Mapping table (mutations passing each filter)

| Cohort | Total SNP | CDS found | Translation match | UniProt found | Window ±5 match | In analysis |
|---|---|---|---|---|---|---|
| George | 34583 | 34408 (99.5%) | 33239 (96.6%) | 31938 (92.4%) | TBD | TBD |
| JHU | 7858 | 0 (RefSeq) | n/a | 6900 (99.0%) | 6051 (87.7%) | TBD |
| CLCGP | 8078 | 8066 (99.8%) | 7784 (96.5%) | 7293 (90.3%) | 6682 (91.9%) | TBD |

### Discard fractions (UniProt + match)
- George: 13.4% (no UniProt 7.6% + mismatch 5.8%)
- JHU: 13.2% (no UniProt 1.0% + mismatch 12.2%)
- CLCGP: 17.3% (no UniProt 9.7% + mismatch 7.6%)

### JHU discrepancy (7858 vs 6972)
- 7858 = all coding SNPs from JHU API.
- 6972 = after filtering to those with transcript_id (NM_) and protein position.
- Difference 886 = mutations without transcript or protein position.

### Mismatch reasons (George/CLCGP)
- To be classified: complementary strand, codon phase shift, transcript mismatch.

### Duplicate samples (removed)
- George vs CLCGP: 23 samples (overlap ≥0.5) → keep George, remove CLCGP.
- JHU vs CLCGP: 1 sample (2334192, overlap 0.675) → keep CLCGP, remove JHU.
- After removal: George 107, JHU 41, CLCGP 6. Total 154 (was 178).

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
