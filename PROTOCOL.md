# Paired ascertainment of diabetes eye care in BRFSS

## Question

How much self-reported retinal-photography contact is invisible to a dilation-only eye-care metric, and does the measurement choice change apparent geographic and socioeconomic inequalities? This is a measurement and service-contact study, not a validation of screening quality, an AI implementation study, or a treatment-effect analysis. It is distinct from the earlier receptor atlas and incretin eye-examination projects.

## Data and eligibility

Use 2025 BRFSS as the primary cross-sectional sample and 2022 as a separate earlier-year replication. Include adults with diagnosed nongestational diabetes in US states/DC administering the Diabetes module on their full combined questionnaire, as explicitly documented by CDC in the annual module-by-state table. Exclude territories and states that require a version-specific module file. This restriction is set for a transparent, reproducible analytic population; estimates represent the participating states, not all US adults. Retain all respondents from eligible states for domain variance estimation. Apply the final full-sample _LLCPWT, design stratum _STSTR and primary sampling unit _PSU, with state identifiers incorporated in design identifiers. Do not pool years. Describe common-state results separately to avoid attributing composition changes to time.

## Fixed definitions and estimands

Annual dilation and annual retinal photography each mean response within the past month or within the past year to the corresponding recency item. Complete pairs classify respondents as both, dilation only, photography only, or neither. The primary estimand is the absolute incremental contact prevalence from adding photography to dilation (the photography-only fraction of all complete-pair respondents). Also estimate the fraction of any reported eye contact missed by dilation alone, and the fraction of respondents classified as lacking dilation who report photography. These are survey-report classifications; photography does not establish retinal image quality, interpretation, referral completion, nonmydriatic technique, AI use, or guideline-concordant screening.

Repeat at a two-year threshold. Do not label nonreceipt as missed indicated screening because disease findings and recommended intervals are unmeasured. Report the conventional dilation prevalence and the composite prevalence side by side with paired covariance. No P value is needed for the mathematical fact that adding a modality cannot reduce prevalence.

## Secondary comparisons and models

Prespecified groups: metropolitan/nonmetropolitan residence; age 18–44, 45–64 and ≥65; sex; race/ethnicity (non-Hispanic White, non-Hispanic Black, Hispanic, other/multiracial); education (<high school, high school, some college, college graduate); annual household income (<$35,000, $35,000–74,999, ≥$75,000, unknown); cost-related inability to see a doctor; serious vision difficulty; current insulin use; reported diabetes type. Describe each four-way pattern in all groups, retaining nulls.

Estimate adjusted associations with photography-only contact using survey-weighted logistic regression followed by marginal standardization over the observed complete-case covariate distribution. Primary inequality contrasts are nonmetropolitan versus metropolitan and cost barrier versus no cost barrier under both dilation-only and either-modality definitions. Adjust for age group, sex, race/ethnicity, education, insulin use and state. Income and serious vision difficulty enter an extended model. Differences between definitions use within-respondent covariance, not independent-sample tests. The four primary adjusted inequality contrasts (two exposures by two definitions) receive Holm correction; all other subgroup results are descriptive. Do not infer causal effects of residence, affordability, imaging or insulin.

## Missingness and precision

Primary complete-pair analysis; report unknown/refused/missing separately for both items. Sensitivities: one positive modality sufficient even if the other is unknown; lower and upper assignments for incremental classification among missing pairs; two-year window; reported type 2 diabetes only; leave-one-state-out; common full-module states across years. No outcome-driven covariate selection. Continuous survey weights are not integer patient counts. Use Taylor linearization of ratio estimands and weighted estimating equations, PSU totals centered within design strata, finite-sample stratum factors and t reference with design degrees of freedom. For singleton strata use a documented conservative grand-mean centering sensitivity; disclose their number and variance handling. Suppress inference where unweighted group n<50 or event count<10; preserve descriptive counts and machine-readable estimates with reliability flags.

## Outputs and interpretation

Retain raw archives, source documentation, SHA-256 manifest, variable coding ledger, all cross-tabulations, full model coefficients/covariance, diagnostics and reproducible code. Produce an English manuscript under 5000 main-text words, tables, high-resolution 600-dpi and vector scientific figures, and a supplement. All results are reported regardless of statistical significance. Any later deviation or data-availability change must be logged. Search novelty against paired-modality BRFSS analyses and related screening-surveillance literature; do not guarantee global priority or editorial acceptance.
