# Prior exploratory analysis

This file records what was analysed before the pre-analysis plan ([`PREREGISTRATION.md`](../PREREGISTRATION.md), tagged `v0.1-prereg`) was written. It expands Section 2 of the plan.

## The earlier draft

A descriptive draft, "Does heat leave a mark? Attention decay and residue after UK heatwaves, 2022–2026" (dated 5 November 2026), was written before the plan. It was not posted. Its analyses are superseded by this repository.

## Data the draft used

| Source | Series |
|---|---|
| Google Trends, UK (national) only, daily and weekly | "heatwave", "climate change", "air conditioning", "heat pump", "net zero" |
| Wikipedia, global English pageviews | Climate change, Heat wave, Air conditioning, Heat pump, Net-zero emissions |
| Heatwave dates | Four episodes dated by hand: 2022, 2023, 2025, 2026 |

## Methods the draft used

- A 60-day pre-onset baseline for each episode.
- A 7-day centred moving average.
- An exponential decay fit with an asymptote, giving a half-life for each episode and series.
- Moving-block bootstrap confidence intervals.
- Comparison of attention 61–90 days after onset with the baseline, with and without seasonal adjustment.

## What the draft found

1. Attention to heat decayed quickly, with half-lives of a few days.
2. "Climate change" searches looked elevated 61–90 days later without seasonal adjustment, but not after adjusting for season.
3. Attention in 2022 did not clearly differ from June 2026.

## Known problems with the draft, addressed by the plan

- **Incomplete 2026 dating.** The draft counted heatwaves in May and 16–30 June 2026 only. Heatwave conditions also occurred in July and August 2026, so its 2026 decay and residue windows overlapped later heat. The plan dates heatwaves by rule from temperature data instead (Section 4).
- **No counterfactual.** The draft compared each episode with its own pre-onset baseline. The plan compares with the same calendar weeks in other years and adds placebo dates (Sections 7–8).
- **No regional variation.** The draft used national series only. The plan adds a four-nation dose-response analysis (Section 7.2).

## Data that were new when the plan was written

Listed in Section 2 of the plan: Trends series for the four UK nations; the terms "global warming", "portable air conditioner" and "heat stroke"; the placebo and anchor terms; all repeated pulls; UK-only Wikipedia pageviews; the Wikipedia articles Met Office, Climate Change Act 2008, Climate Change Committee, Heat stroke, Fan (machine) and Chess; and all temperature data.
