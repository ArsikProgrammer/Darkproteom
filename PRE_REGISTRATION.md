# Pre-registration: Analysis Plan (Dark Proteome × Tobacco Mutational Signature)

**Status:** Fixed before viewing real results.  
**Date:** 2026-10-25  
**Project:** Dark proteome and mutational footprint of tobacco smoke

---

## 1. Hypothesis

**H1 (working).** In tumors with stronger smoking mutational signature (SBS4), the fraction of mutations falling in intrinsically disordered regions (IDR) of proteins is higher than expected given the tumor's own 96-context mutation spectrum and the distribution of those contexts between dark and ordered protein positions. The effect is expected to be stronger in SCLC than in LUAD/LUSC.

**H0.** After correction for context composition, coverage, and selection, no systematic shift remains: observed/expected ratio of mutations in IDR is consistent with 1.

---

## 2. Primary definition of "darkness"

- **Primary:** AlphaFold pLDDT < 50 (very low confidence / likely disordered).
- **Sensitivity:** pLDDT < 30 and pLDDT < 70.
- The primary threshold is fixed before analysis; sensitivity thresholds are exploratory only.

---

## 3. Mutation-to-residue mapping (blocking step)

- Transcript: GENCODE v19, protein-coding CDS.
- Translate CDS, align to UniProt sequence.
- Keep only mutations where translated residue matches UniProt residue.
- Record and report the fraction of mutations discarded due to mismatch.
- Validate on known mutations (TP53 R175H, KRAS G12C, and others from our data where reference amino acid must match).
- Check reference allele matches hg19; amino acid matches UniProt.
- Do not silently fix mismatches; count them.
- **Do not examine is_dark vs SBS4 during this step.**

---

## 4. Samples and cohorts

- Primary: SCLC (George, JHU, CLCGP).
- Secondary: TCGA LUAD and LUSC (dose-response, non-smokers).
- Check sample overlap between George, CLCGP, JHU; remove duplicates.
- Minimum mutations per sample for signature estimation: determined by bootstrap stability.

---

## 5. Signatures

- Primary predictor: **SBS4 fraction** (SBS4 / total mutations), standardized (1 SD).
- Fitting set: SBS1, SBS2, SBS4, SBS5, SBS13, SBS17, SBS40.
- WGS (George): fit on whole genome.
- WES (JHU, CLCGP, TCGA): fit on exome; normalize to exome.
- Report bootstrap CI for SBS4.
- SBS40 is exploratory only (COSMIC: etiology unknown).

---

## 6. Primary model (fixed)

**Level: mutation (binomial).**

logit P(is_dark=1) = β0 + β1 · SBS4_std + β2 · TMB_std + (1 | dataset)

- **Main coefficient:** β1 — effect of 1 SD of SBS4 fraction on log-odds of a mutation being in IDR.
- **Primary effect measure:** OR per 1 SD SBS4 = exp(β1).
- Dataset as random effect.
- Synonymous mutations as negative control.
- SBS1 (age) as negative control: should not predict IDR shift.

---

## 7. Equivalence boundary (fixed before analysis)

- **Minimum interesting effect:** OR = 1.2 per 1 SD SBS4.
- **Equivalence interval:** (1/1.2, 1.2) = (0.833, 1.2).

**Interpretation rules:**
- If 95% CI for OR lies entirely within (0.833, 1.2) → **informative null** (effect smaller than of interest).
- If 95% CI is wider and crosses boundary → **inconclusive** (cannot exclude meaningful effect).
- If CI excludes 1 and lies outside interval → **effect detected**.

This prevents calling every narrow CI "inconclusive" and every wide CI "null".

---

## 8. O/E calculation

- Observed: mutations in IDR.
- Expected: based on tumor's own 96-context spectrum and position-level context/darkness map.
- Metric: log(O/E).
- Absolute position counts (not gene fractions).

---

## 9. Covariates (after primary analysis)

- Gene expression.
- Replication timing.
- (Chromatin deferred.)

---

## 10. Exclusion criteria (fixed)

- Samples with fewer than N mutations (N set by bootstrap stability).
- Duplicate samples across cohorts (remove, record count).
- Mutations not mappable to protein residue (record discard fraction).
- If SBS4 is unstable (bootstrap CI includes 0 or very wide), record as limitation and proceed with caveat.

---

## 11. Sensitivity analyses

- Thresholds: pLDDT < 30, < 50, < 70.
- Signature set variations.
- Per-cohort analysis.

---

## 12. Interpretation rules

- No causal claims.
- If SBS4 CI includes zero but is narrow → informative null (per equivalence boundary).
- If SBS4 CI is wide → inconclusive.
- Report power analysis before collecting more data.

---

## 13. Synthetic data policy

- All synthetic/test files removed from analysis directory.
- No synthetic data used in final results.

---

## 14. Commit

- This file committed to git with date before any real analysis.
- Optionally uploaded to OSF/GitHub release for public timestamp.
