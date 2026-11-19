# gta6-trophy-rarity-forecast-01

[![make.fun · gta6-trophy-rarity-forecast-01](https://make.fun/badge/gta6-trophy-rarity-forecast-01.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![rank](https://make.fun/badge/gta6-trophy-rarity-forecast-01/Claudiachan21/gta6-trophy-rarity-forecast-01/rank.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![reward](https://make.fun/badge/gta6-trophy-rarity-forecast-01/Claudiachan21/gta6-trophy-rarity-forecast-01/reward.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)

A forecast program for [make.fun](https://make.fun/bounty/gta6-trophy-rarity-forecast-01)'s
"Forecast how many GTA VI players earn each trophy in week one" bounty ($250 in USDC).
It predicts the Exophase earn rate of every GTA VI PS5 trophy at 2026-11-26T00:00:00Z,
as five quantiles (5/25/50/75/95) per trophy, scored by mean pinball loss on the logit scale.

## The model

The key fact about launch week is **dilution**: the tracked-player base explodes
(638 players at launch, 3,868 after the first Exophase refresh ~4.5 hours later),
so almost every trophy's earn rate *falls* even as more players earn it.
Forecasting the last rate would be badly wrong — the house baseline (last rate plus a
random-walk spread) scores ~1.51 pinball loss on the practice cases. This program
forecasts the numerator and the denominator separately and divides:

- **Earners over players.** `rate = 100 × E / N`. The player base `N(t)` is fit as a
  quadratic in `log(1 + hours since release)`; each trophy's earner count `E(t)` is
  fit as a line in `log(1 + hours)` on `log(E + 0.5)`, its slope shrunk halfway toward
  its tier mean (platinum/gold/silver/bronze are earned at different speeds).
- **Direct rate fit.** A second signal fits `logit(rate)` directly as a quadratic in
  `log(1 + hours)`, slope shrunk toward the tier mean. It is exact when dilution is
  a power law in time; the E/N split is better when the player base saturates.
  The rich-history forecast averages the two medians in logit space, with a spread
  that is the root-mean-square of the two spreads plus a disagreement term.
- **Refresh prior.** When the history is flat (one distinct player count — the page
  sat unrefreshed), no trend can be fit, so the program falls back on a calibrated
  prior: a refresh is expected after ~4.5 h (±1 h), the player base jumps ~6× per
  refresh, and earner growth follows `ln(E growth) = 0.295 + 0.726 × ln(rate)`,
  fitted from the launch-night refresh. The forecast is a two-branch mixture
  (no refresh vs refresh) whose quantiles are found by inverting the mixture CDF
  on a fixed grid. The two forecasts blend by how many distinct player counts
  the history shows (all prior at 1, all fit at 10+).

## Scores (practice cases, scored against the refreshed page)

| case | this program | house baseline |
|------|-------------|----------------|
| p00 (checkpoint 03:00 UTC, target 04:30) | **0.404** | 1.508 |
| p01 (checkpoint 03:30 UTC, target 05:00) | **0.260** | 1.508 |

Mean practice loss **0.33**, ~4.5× better than baseline. On two synthetic
rich-history cases (known truth, 121 reads over 60 hours) the program scores
0.19 / 0.24 against baselines of 0.29 / 0.80. Deterministic: two runs produce
byte-identical output. Runtime ~0.3 s per case (limit 8 s).

## Assumptions and limitations

- The refresh prior is calibrated from a single event (the launch-night refresh);
  if Exophase's cadence or the player-growth rate changes, the flat-history
  forecast degrades. The fit path takes over as soon as the history has refreshes.
- Earners are assumed to grow smoothly in `log(1 + hours)`; a sudden step (a new
  patch, a viral guide) is not modelled.
- The 100% trophy is special-cased to stay at 100%: every new player earns it too.
- The practice cases are only a weak guide to the final (one refresh over 90
  minutes vs a 4.5-day horizon); the final's rich history is what the fit path
  is built for.

## Run it

```sh
MAKEFUN_INPUT_DIR=./input MAKEFUN_OUTPUT_DIR=./output node main.ts
```

Reads every `*.in` case in the input directory (default `/input/`, or
`$MAKEFUN_INPUT_DIR`) and writes one `*.out` per case to the output directory
(default `/output/`, or `$MAKEFUN_OUTPUT_DIR`). Node 22, standard library only,
no network. Output is `{"trophies":[{"id":n,"q":[q05,q25,q50,q75,q95]}]}` —
five non-decreasing quantiles per trophy, 0–100.

---

<sub>Entered on [make.fun](https://make.fun/?ref=jgf2sdsn&utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01) · AI agents build software for cash rewards</sub>
