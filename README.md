# Does heat leave a mark?

**Causal effects of UK heatwaves on public attention, and their decay and residue, 2022–2026**

Kaustubh Ramekar · independent researcher, London

## Question

When a heatwave hits the UK, attention to heat and climate change rises. This project asks the causal version of that observation:

- How much does heatwave exposure raise attention, compared with what the same weeks would have shown without the heat?
- How quickly does the effect decay, and is anything left two to three months later (*residue*), or does attention fully *reset*?
- Do *coping* searches (staying cool) and *concern* searches (climate change, net zero) behave differently?
- Within the same week, do hotter UK nations respond more?
- Has the response weakened since the record-breaking summer of 2022?

## Design

- **Treatment:** heatwave days defined by rule from HadUK-Grid daily maximum temperature (ERA5 for 2026), following the Met Office three-day definition.
- **Outcomes:** Google Trends (UK and the four nations, anchored and averaged over repeated pulls) and Wikipedia pageviews (UK-only and global).
- **Analysis 1:** national event study across five summers, with week-of-year and year fixed effects, onset- and end-aligned.
- **Analysis 2:** four-nation dose-response with week fixed effects, so identification comes only from differences in heat between nations in the same week.
- **Checks:** pre-trend leads, placebo dates, placebo years, placebo outcomes, anchor validity and concurrent-event exclusions.

```
season ──► heat ──► attention ◄── year
            │           ▲
            └─► news ───┘      (mediator, not controlled)
```

## Pre-registration

The full analysis plan is in [`PREREGISTRATION.md`](PREREGISTRATION.md). It was the first commit of this repository and is tagged [`v0.1-prereg`](../../releases/tag/v0.1-prereg), before any outcome data were downloaded. Any later change is logged in [`docs/DEVIATIONS.md`](docs/DEVIATIONS.md), each committed before the analysis it affects was run:

1. **Pull dates.** Google Trends pulls 2 and 3 ran on 9 and 10 November instead of 7 and 8 November; the three-different-days rule still holds.
2. **Sparsity rule.** As written, the rule removes "heatwave" and every coping term, because these searches are near zero outside summer. The plan-as-written results are reported, alongside a labelled analysis that keeps all pre-registered terms, with a safeguard against floor effects.
3. **Calibrated inference.** A size check on simulated noise shows the pre-registered Newey–West and Driscoll–Kraay tests reject far more often than 5% with so few heatwaves. Tests using the placebo and permutation distributions are added; the pre-registered tests are still reported.

What was analysed before the plan was written is described in [`docs/PRIOR_EXPLORATORY_ANALYSIS.md`](docs/PRIOR_EXPLORATORY_ANALYSIS.md).

## Status

| Step | Status |
|---|---|
| Pre-analysis plan | Done (5 Nov 2026) |
| Configuration (event dates, concurrent events) | Done (5 Nov 2026) |
| Temperature data and heatwave rule | 2021–2025 done (HadUK-Grid); 2026 from ERA5, calibration years downloading |
| Google Trends pulls (3 days) | Pull 1 done (6 Nov, 15:13–15:58 UTC); pull 2 done (9 Nov, 14:38–14:46 UTC); pull 3 planned 10 Nov (see deviation 1) |
| Wikipedia pageviews | Done (UK-only and global) |
| Analyses 1 and 2, placebos | Code done (`src/06`–`src/09`); final run after pull 3 |
| Paper | Not started |

## Repository layout

```
PREREGISTRATION.md   analysis plan (fixed before outcome data)
config/              event dates and concurrent events, fixed in advance
data/raw/            raw downloads, never edited
data/processed/      cleaned panels built by src/
data/MANIFEST.csv    source, parameters, download time and SHA-256 of each raw file
docs/                deviations log and prior exploratory analysis
src/                 download, cleaning and analysis scripts
figures/ results/    outputs, all produced by `make all`
tools/               helper scripts
```

## Reproducing

```
conda env create -f environment.yml
conda activate heat-leaves-a-mark
make all
```

## Data sources

- Met Office HadUK-Grid (via CEDA), Open Government Licence
- Copernicus ERA5 post-processed daily statistics (Climate Data Store), CC-BY 4.0. Contains modified Copernicus Climate Change Service information.
- Natural Earth 1:10m map units, public domain
- ONS mid-2024 population estimates, Open Government Licence
- Google Trends
- Wikimedia pageviews (REST API and differentially private country data), CC0
