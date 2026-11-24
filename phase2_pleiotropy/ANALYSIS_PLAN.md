---

editor_options: 
  markdown: 
    wrap: 72
---

# Phase 2 analysis plan: fibromyalgia, genetic overlap, and causal direction

**Status:** DRAFT v0.1. To be finalized and committed (and git-tagged `plan-v1`) **before any Mendelian randomization (MR) is run**. Any change made after seeing MR results will be listed under "Deviations" at the bottom, with the reason.

## Question

Fibromyalgia's genetics overlap heavily with pain, psychiatric and metabolic traits. The source GWAS (Kerrebijn et al. 2026) raises, but does not resolve, whether that overlap reflects a shared central nervous system vulnerability or effects that arrive secondarily through comorbidities such as obesity.

This phase asks, in three parts:

1.  **Genetic correlation.** Do we reproduce the paper's correlations with pain, psychiatric and digestive traits?

2.  **Causal direction (MR).** Does genetic liability to a candidate exposure predict fibromyalgia risk?

3.  **Independence (multivariable MR).** Which exposures retain an effect when considered jointly?

## Data and overlap handling

- **Fibromyalgia GWAS:** the paper's European-only analyses (GWAS Catalog GCST90838603-GCST90838628). The study provides leave-one-cohort-out versions. For every exposure we use the version that removes the cohort(s) it shares people with:
  - FinnGen-derived traits: `GCST90838609` (European, leave-FinnGen-out).
  - UK-Biobank-derived traits: `GCST90838615` (European, leave-UK-Biobank-out).
  - Other exposures: leave-out version chosen after checking each study's cohort list (to confirm).
- **Exposure / correlation traits** (sources to be verified for open access, ancestry and build): low back pain, IBS, depression, insomnia, PTSD (FinnGen R13 endpoints for genetic correlation); BMI, multisite chronic pain, insomnia, depression (large GWASs for MR).
- Residual overlap that cannot be removed will be stated explicitly in the write-up.

## Part 1: genetic correlation

- LD score regression through GenomicSEM (R), European LD scores, HapMap3 SNPs.
- Validation target: reproduce the paper's heritability (done: 11.6% vs 10.4%) and the direction and rough size of its correlations (e.g. low back pain 0.75, PTSD 0.78, IBS 0.70, depression 0.63, insomnia 0.47).
- Multiple testing: Bonferroni across all correlations tested in this part.

## Part 2: Mendelian randomization

- **Direction:** exposure -\> fibromyalgia is primary. Fibromyalgia -\> exposure is secondary and only attempted where instruments are adequate.
- **Instruments:** P \< 5e-8, LD-clumped (r2 \< 0.001, 10,000 kb window, European reference), single nucleotide variants; palindromic variants with ambiguous strand dropped. Instrument strength reported as per-SNP and overall F-statistics; weak instruments (F \< 10) noted.
- **Primary estimator:** inverse-variance weighted (random effects).
- **Sensitivity estimators:** MR-Egger, weighted median, weighted mode; Cochran's Q heterogeneity; MR-Egger intercept; leave-one-variant-out; Steiger directionality.
- **Multiple testing:** Bonferroni across the number of exposures tested (k to be fixed here).
- **Decision rule (fixed in advance):**
  - *Supportive evidence* = primary IVW P \< 0.05 / k, same direction in at least two of the three sensitivity estimators, and no significant MR-Egger intercept.
  - *Inconclusive* = anything else. Inconclusive results are reported in full, not dropped.

## Part 3: multivariable MR

- Candidate set: BMI, insomnia, depression, multisite chronic pain (final set depends on instrument availability after Part 2).
- Conditional F-statistics reported. Interpretation limited to which exposures keep an independent effect, not to magnitudes of causal effect.

## Reporting commitments

- All tested exposures are reported, including null results.
- Forest plots and a table of every estimate, with the number of instruments used.
- Limitations stated: sample overlap, instrument pleiotropy, ICD-coded fibromyalgia, European-only ancestry, and the usual MR assumptions.
- No claim of causation: results are described as evidence consistent or inconsistent with a causal effect under the assumptions.

## To confirm before freezing

- [ ] Source, ancestry, genome build and cohort list for each MR exposure GWAS
- [ ] Which fibromyalgia leave-out file pairs with each exposure
- [ ] How clumping will be done (local reference panel vs an online service)
- [ ] The value of k (number of exposures) for the Bonferroni threshold

## Deviations

*(none yet)*
