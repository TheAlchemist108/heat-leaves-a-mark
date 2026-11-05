# Pre-analysis plan

**Study:** Does heat leave a mark? Causal effects of UK heatwaves on public attention, and their decay and residue, 2022–2026
**Author:** Kaustubh Ramekar (independent researcher, London, UK)
**Plan version:** 1.0, written 5 November 2026
**Status when written:** No data for the analyses in this plan have been downloaded. Some overlapping data were analysed earlier for a descriptive draft; Section 2 says exactly what was seen.
**How this plan is time-stamped:** This file is the first commit of the public repository and is tagged `v0.1-prereg` with a GitHub release. An independent snapshot is requested from Software Heritage (`tools/swh_save.sh`). The SHA-256 of this file is recorded in the release notes. Later changes are made only as new commits and are logged in `docs/DEVIATIONS.md`. History is never rewritten.

---

## 1. Questions

Earlier descriptive work asked whether attention rises and falls after UK heatwaves. This plan asks the causal version of those questions:

- **Q1 (effect and decay).** How much does heatwave exposure raise public attention, measured against what the same calendar weeks would have shown without the heat, and how fast does the effect decay?
- **Q2 (residue).** Is anything left 9–12 weeks after the heat (*residue*), or does attention return to its counterfactual level (*reset*)?
- **Q3 (mechanism).** Do *coping* searches (staying cool) and *concern* searches (climate change and policy) behave differently?
- **Q4 (dose).** Within the same week, do UK nations that got hotter show larger responses?
- **Q5 (habituation).** Per unit of heat, is the response smaller in 2023–2026 than in 2022?

The estimand throughout is the **total effect** of heat exposure on attention, including effects that work through news coverage. News coverage is a mediator and is not controlled for (Section 5, DAG).

## 2. What was known before this plan (disclosure)

This plan is **not** written blind to all data. A descriptive draft ("Does heat leave a mark? Attention decay and residue after UK heatwaves, 2022–2026", dated 5 November 2026, not posted) analysed:

- national (UK) Google Trends daily and weekly series for "heatwave", "climate change", "air conditioning", "heat pump" and "net zero";
- global English Wikipedia pageviews for "Climate change", "Heat wave", "Air conditioning", "Heat pump" and "Net-zero emissions";
- four hand-dated episodes (2022, 2023, 2025, 2026).

Its main findings were: attention to heat decayed with half-lives of a few days; "climate change" searches looked elevated 61–90 days later without seasonal adjustment, but not after adjusting for season; and 2022 did not clearly differ from June 2026.

Two consequences for this plan:

1. **Predictions below are informed by that draft.** Where the draft already gives a likely answer (Q2 for "climate change"), the prediction is stated, and the test is designed so that a null result is informative (an equivalence test, Section 7.3).
2. **The draft's 2026 dating was incomplete.** It counted heatwaves in May and 16–30 June 2026 only. Heatwave conditions also occurred in July and August 2026, so its S2026 decay and residue windows overlap later heat. This plan replaces hand-dated episodes with a rule applied to temperature data (Section 4.1).

**Genuinely new, unseen data in this plan:** Trends series for the four UK nations; the new terms "global warming", "portable air conditioner", "heat stroke", the placebo terms and both anchor terms; all repeated pulls; UK-only Wikipedia pageviews (Wikimedia differentially private country data); the Wikipedia articles "Met Office", "Climate Change Act 2008", "Climate Change Committee", "Heat stroke", "Fan (machine)" and "Chess"; and all temperature data.

Full details are in `docs/PRIOR_EXPLORATORY_ANALYSIS.md`.

## 3. Study period and units

- **Primary sample:** weeks from Sunday 7 November 2021 to Saturday 31 October 2026 (260 weeks). This covers five summers: 2022, 2023, 2024, 2025 and 2026.
- **Weeks** run Sunday to Saturday, matching Google Trends weekly data.
- **Units:** the UK (national analyses) and England, Scotland, Wales and Northern Ireland (dose-response).
- **Daily windows** (for descriptive decay curves only): national, 1 April to 15 December each year; 1 April to 3 November for 2026.

## 4. Treatment: heat exposure

### 4.1 Temperature data

- **2021–2025:** HadUK-Grid v1.3.2.ceda, daily maximum air temperature, official area averages for each UK country (`country/tasmax/day`, 1931–2025, Met Office via CEDA, Open Government Licence).
- **2026:** HadUK-Grid does not yet cover 2026. Daily maximum 2 m temperature comes from ERA5 (Copernicus Climate Data Store, "ERA5 post-processed daily statistics on single levels"), averaged over each nation's land area using Natural Earth 1:10m map-unit boundaries.
- **Calibration:** For each nation, ERA5 nation averages are regressed on HadUK-Grid nation averages for April–October 2016–2025 (daily, linear). 2026 values are the calibrated ERA5 values.
  - **Rule:** if the daily correlation is below 0.95 for a nation, 2026 is excluded from that nation's dose-response analysis, and from the national analysis if the population-weighted UK series fails.
- **UK series:** a population-weighted mean of the four nation series. The weights come from the ONS mid-2024 population estimates and are held fixed for all years.

### 4.2 Heatwave days and onsets (primary, by rule)

For each series s (each nation, and the UK):

- **Threshold:** T*_s is the 95th percentile of daily maximum temperature over June–August 1991–2020.
- **Heatwave day:** a day belonging to a run of at least 3 consecutive days with Tmax ≥ T*_s. The three-day rule follows the Met Office heatwave definition. The percentile threshold replaces the Met Office's county thresholds, which are not published as a gridded product.
- **Heatwave spell:** consecutive heatwave days. Spells separated by fewer than 14 days merge into one **episode**.
- **Onset week:** the week containing the first day of an episode.
- **Heat dose (nation n, week w):** the sum over the week's days of max(0, Tmax − T*_n), in degree-days.

### 4.3 Alternative treatment definitions (robustness)

- Thresholds at the 90th and 97.5th percentiles.
- Met Office–documented heatwave dates (`config/events_metoffice.yaml`).
  - Only dates stated in Met Office event reports, monthly weather reports or climate summaries are used. Where the Met Office gives no dates for an event, that event is omitted from this list rather than dated from other sources.
  - This list is finalised and committed before any outcome data are downloaded.

## 5. Causal assumptions (DAG)

```
season (week of year) ──► heat ──► attention
season ─────────────────────────► attention
year (trends, platform change) ──► attention
heat ──► news coverage ──► attention        (mediator: not controlled)
heat in Europe ──► news ──► attention       (affects the national analysis; absorbed by week fixed effects in the dose-response)
concurrent events (COP, elections, deaths of public figures) ──► attention
Trends sampling, anchor term ──► measured attention
```

**Identifying assumption, national model:** after removing week-of-year and year effects, the timing of heatwave onsets is unrelated to other shocks to attention. Leads (pre-onset coefficients) and placebo dates test this.

**Identifying assumption, dose-response:** within a given week, differences in heat dose between nations are unrelated to differences in other attention shocks, after nation-by-year and nation-by-week-of-year effects. UK-wide news, including news about heat elsewhere, is absorbed by the week fixed effects.

## 6. Outcomes

### 6.1 Google Trends

**Settings:** web search, all categories. Geos: GB (UK), GB-ENG, GB-SCT, GB-WLS and GB-NIR.

**Anchor terms.** Every request contains the anchor term **"dentist"**, so all values can be expressed relative to it. The secondary anchor is **"plumber"**, used only for the anchor validity check (Section 8).

**Term families:**

| Family | Terms | Role |
|---|---|---|
| Salience | heatwave | First stage: was the heat noticed? |
| Coping | air conditioning; portable air conditioner; heat stroke | Coping with heat |
| Concern | climate change; global warming; net zero | Climate concern and policy |
| Transition (exploratory) | heat pump | Not confirmatory |
| Placebo outcomes | car insurance; recipes | Should not respond to heat |
| Anchors | dentist; plumber | Normalisation |

**Requests.** Each request has five terms or fewer:

- **R1:** dentist, heatwave, climate change, air conditioning, car insurance
- **R2:** dentist, global warming, net zero, portable air conditioner, heat stroke
- **R3:** dentist, heat pump, recipes, plumber, climate change. Weekly only; the repeated "climate change" checks consistency across requests.

**Windows:**

- **Weekly:** 2021-11-07 to 2026-10-31, all five geos, R1–R3.
- **Daily:** national (GB) only, R1–R2, using the Section 3 windows for each year.

**Repeated pulls.** Three complete pulls on three different days, at least 20 hours apart. Planned for 6, 7 and 8 November 2026.

**Measurement:**

- Values reported as "<1" are coded 0.5. The log ratio for term j is r_j = ln(max(x_j, 0.5)) − ln(max(x_anchor, 0.5)).
- The outcome is the mean of r_j over the three pulls. The SD across pulls and the between-pull correlation are reported.
- Family indices are the mean of the member terms' log ratios, so effects are in log points.

**Sparsity rule.** This is applied before any effect is estimated. A term-by-geo series is dropped if more than 20% of its weeks have a pull-averaged raw value of 1 or less. Family indices are then built from the remaining terms. If the anchor itself fails in a geo, that geo is dropped.

### 6.2 Wikipedia

**Articles,** identified by page ID:

| Family | Article (page ID) |
|---|---|
| Salience | Heat wave (293270); Met Office (319671) |
| Coping | Air conditioning (7221088); Heat stroke (24999087); Fan (machine) (12877572) |
| Concern | Climate change (5042951); Climate Change Act 2008 (10031054); Climate Change Committee (10031822) |
| Transition (exploratory) | Heat pump (68316) |
| Placebo | Chess (5134) |

**Sources:**

- **(a) UK-only (primary for Wikipedia):** Wikimedia differentially private daily pageviews by country, project and page. Rows with country = GB and project = en.wikipedia are kept.
  - The source has two releases: 9 Feb 2017–5 Feb 2023 (historical) and 6 Feb 2023 onwards (current). The two use different algorithms, and from 25 May 2023 bot filtering changed.
  - Models using this source include step indicators for 6 Feb 2023 and 25 May 2023.
  - Days below the release threshold are missing. An article is dropped if more than 20% of its days are missing.
- **(b) Global (secondary):** Wikimedia REST API per-article daily pageviews (all-access, user agent).
  - Title moves during the sample are checked with the MediaWiki move log and recorded.

**Outcome:** the log of weekly views. Family indices are means of the member articles' log views.

## 7. Analyses

### 7.1 National event study (Analysis 1, primary for Q1–Q3)

For each outcome y (UK, weekly), two alignments are estimated:

```
(1a) onset-aligned:  y_w = Σ_{k=-4}^{12} β_k · Onset_{w−k} + μ_{week-of-year(w)} + τ_{year(w)} + ε_w
(1b) end-aligned:    the same model with End_w (the week containing an episode's last heatwave day) in place of Onset_w
```

**Why two alignments.** Episodes differ in length; the 2026 episode, for example, runs for weeks. Impact is measured from the onset. Decay and residue are measured from the end of the heat, so that "weeks 9–12 later" never falls inside a continuing heatwave. In model (1b) the leads cover the episode itself, so the pre-trend test uses model (1a) only.

- **Specification:**
  - Onset_w = 1 in onset weeks (Section 4.2).
  - Overlapping event windows are handled in distributed-lag form (Schmidheiny and Siegloch, 2023).
  - The window endpoints are binned: k = −4 collects all earlier leads and k = 12 all later lags.
  - The reference period is k = −1.
  - "Year" is the calendar year.
- **Inference:** Newey–West standard errors with 8 lags, plus randomisation inference (Section 8).
- **Estimands:**
  - **Impact:** mean of β_0, β_1, from model (1a).
  - **Decay:** mean of β_4, …, β_8, from model (1b).
  - **Residue:** mean of β_9, …, β_12 (about days 61–90 after the heat ends), from model (1b).
  - **Half-life:** β_k = A·exp(−λk) fitted by nonlinear least squares over k = 0…12 of model (1b); t½ = ln 2/λ. The confidence interval comes from a moving-block bootstrap of the regression residuals (block of 8 weeks, 2,000 replicates). The half-life is reported as "not identifiable" if λ hits its bound or the upper confidence limit exceeds 12 weeks.

**Descriptive daily decay (continuity with the draft).** For each episode and each daily national series, the draft's method is repeated: a 60-day pre-onset baseline, a 7-day centred mean, an exponential fit with asymptote, and a moving-block bootstrap. Episodes are now dated by rule. Where fewer than 60 days of baseline fall inside the daily window, the available days are used and the result is flagged. These results are descriptive and are not part of the confirmatory tests.

### 7.2 Four-nation dose-response (Analysis 2, primary for Q4–Q5)

For outcome y in nation n and week w (regional Trends; Wikipedia has no regional breakdown):

```
y_{n,w} = Σ_{k=0}^{12} β_k · Dose_{n,w−k} + α_{n,year} + γ_{n,week-of-year} + λ_w + ε_{n,w}
```

- **Fixed effects:** λ_w (week fixed effects) absorbs everything common to the UK in a given week. Identification comes only from differences in heat between nations within the same week.
- **Estimands:**
  - **Short-run:** β_0 + β_1, per 10 degree-days.
  - **Residue:** mean of β_9…β_12.
  - **Habituation:** add Dose × 1[year ≥ 2023] terms for k = 0, 1. The estimand is the interaction on the short-run sum.
- **Inference:**
  - Driscoll–Kraay standard errors (8 lags).
  - Randomisation inference by permuting whole nation dose series across nations. With 4 nations there are 24 permutations, so the smallest attainable p-value is 1/24 ≈ 0.042.
  - Year-shift placebo: each nation's dose is moved to the same weeks of another year.
- **Family composition:** in this analysis, a family index uses only the terms that pass the sparsity rule in all four nations, so the index means the same thing in each nation.
- **Limits stated in advance:** 4 nations is a small number, and England dominates. This analysis corroborates the national timing result and is not the headline.

### 7.3 Confirmatory hypotheses

Seven tests, reported with Holm adjustment as one family.

| # | Hypothesis | Outcome | Estimand | Prediction |
|---|---|---|---|---|
| H1 | Heat raises salience | Trends "heatwave" (UK) | Impact | Positive |
| H2 | Coping attention decays within about 2 months | Coping index (UK) | Decay | Equivalent to zero |
| H3 | Concern attention leaves no residue | Concern index (UK) | Residue | Equivalent to zero (from the draft) |
| H4 | Coping responds more than concern | Coping minus concern (UK) | Impact difference | Positive |
| H5a | Hotter nations respond more (coping) | Coping index (4 nations) | Short-run | Positive |
| H5b | Hotter nations respond more (concern) | Concern index (4 nations) | Short-run | Two-sided, no prediction |
| H6 | Habituation after 2022 | Coping index (4 nations) | Interaction | Two-sided; the draft suggests little difference |

**Decision rules:**

- **Holm adjustment** is applied to the model-based p-values: Newey–West (Analysis 1) or Driscoll–Kraay (Analysis 2). For H2 and H3, the p-value entering Holm is the equivalence-test p-value, which is the larger of the two one-sided p-values.
- **Effect present:** the Holm-adjusted p-value is below 0.05 *and* the unadjusted randomisation-inference p-value is below 0.05.
- **Equivalence (H2, H3):** two one-sided tests at α = 0.05 with bounds of ±0.10 log points (about ±10%). Residue is labelled **reset** only if equivalence holds. It is labelled **residue** if the effect is positive by the rule above, **undershoot** if negative by that rule, and **inconclusive** otherwise.
- **Minimum detectable effect:** reported for each test as 2.8 × the SD of the placebo-date distribution, so readers can judge what a null rules out.

Everything else (single terms, Wikipedia, the transition family, daily decay curves, per-episode estimates) is secondary or exploratory. These are reported in full, with unadjusted 95% CIs, and labelled as such.

### 7.4 Habituation (Q5), secondary national version

From Analysis 1, impact is estimated separately for each episode by interacting Onset with episode indicators. Episode impact is then regressed on episode peak dose with a 2022 indicator. With about 5–10 episodes this is descriptive. The 2022 versus 2026 contrast is reported because 2026 was the episode closest in severity.

## 8. Placebo and robustness checks (all pre-specified)

1. **Leads:** a joint Wald test of β_{−4}, β_{−3}, β_{−2} = 0 in model (1a).
2. **Placebo dates:** 2,000 draws of fake onset sets, each with the same number of onsets as the real data. Fake onsets are drawn uniformly from April–September weeks with no real heatwave day within ±4 weeks. The randomisation-inference p-value is the share of draws with |estimate| ≥ |real estimate|.
3. **Placebo years:** real onsets are shifted by ±1 and ±2 years to the same week of year. Shifted onsets within 4 weeks of a real heatwave are dropped.
4. **Placebo outcomes:** "car insurance", "recipes" and Wikipedia "Chess" should show no impact.
5. **Anchor validity:** Analysis 1 is run on ln(dentist) − ln(plumber) from R3. If the impact estimate has p < 0.05, the anchor is treated as heat-sensitive, and all Trends results are re-run with "plumber" as the anchor. Both sets are reported, and the plumber-anchored set becomes primary.
6. **Concurrent events:** re-run with these weeks dropped: COP27, COP28, COP29 and COP30 (from 3 days before the opening to the close); 8–19 September 2022 (death and funeral of Queen Elizabeth II); and the 2024 general election week. Dates are in `config/concurrent_events.yaml`.
7. **Threshold:** the 90th and 97.5th percentiles. **Dating:** Met Office–documented dates.
8. **Pull noise:** results from each single pull compared with the three-pull mean.
9. **Global vs UK Wikipedia:** both sources compared on the overlapping articles.

## 9. Data exclusions and missing data

- **Partial data:** a week is incomplete if Google flags it as partial. Incomplete weeks at the end of the sample are dropped.
- **Wikipedia:** rules as in Section 6.2.
- **Temperature:** missing temperature data are not imputed. If any are missing, the affected weeks are dropped and reported.
- **2026 lags:** for episodes close to the sample end, lags beyond the last observed week are not observed. These episodes add information only to the lags they reach.

## 10. What would change my conclusions

- **Residue in concern attention after heat** would need the H3 rule to show "residue", corroborated by a positive dose-response residue (Analysis 2) and by the UK-only Wikipedia concern index.
- **A failed identification check** (leads, placebo dates or placebo outcomes showing effects of similar size to the real ones) means the design cannot separate heat from coincident shocks for that outcome. The paper will say so, and will not report the estimate as causal.

## 11. Deviations

Every change after this plan is committed is logged in `docs/DEVIATIONS.md`. Each entry records the date, what changed, why, and whether any outcome data had been seen at that point. The paper reports all deviations.

## 12. Software and reproducibility

- Python versions are pinned in `environment.yml`.
- Raw downloads are never edited. Each raw file's SHA-256, source URL or request parameters, and download time are written to `data/MANIFEST.csv`.
- `make all` reproduces every table and figure from the raw files with a fixed random seed (20261105).

## References

- Met Office. Heatwave definition. https://weather.metoffice.gov.uk/learn-about/weather/types-of-weather/temperature/heatwave
- Met Office / CEDA. HadUK-Grid v1.3.2.ceda. https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb
- Copernicus Climate Change Service. ERA5 post-processed daily statistics on single levels. Climate Data Store.
- Wikimedia Foundation. Differentially private pageviews (country–project–page). https://analytics.wikimedia.org/published/datasets/country_project_page/00_README.html
- Schmidheiny, K. and Siegloch, S. (2023). On event studies and distributed-lags in two-way fixed effects models: identification, equivalence, and generalization. *Journal of Applied Econometrics* 38(5), 695–713.
- Driscoll, J. C. and Kraay, A. C. (1998). Consistent covariance matrix estimation with spatially dependent panel data. *Review of Economics and Statistics* 80(4), 549–560.
- Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica* 55(3), 703–708.
- Lakens, D. (2017). Equivalence tests: a practical primer for t tests, correlations, and meta-analyses. *Social Psychological and Personality Science* 8(4), 355–362.
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics* 6(2), 65–70.
