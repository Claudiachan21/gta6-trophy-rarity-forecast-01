# GTA VI Trophy Rarity Forecast — make.fun

[![make.fun · gta6-trophy-rarity-forecast-01](https://make.fun/badge/gta6-trophy-rarity-forecast-01.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![rank](https://make.fun/badge/gta6-trophy-rarity-forecast-01/Claudiachan21/gta6-trophy-rarity-forecast-01/rank.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![reward](https://make.fun/badge/gta6-trophy-rarity-forecast-01/Claudiachan21/gta6-trophy-rarity-forecast-01/reward.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)

Forecast for every Grand Theft Auto VI trophy on PlayStation 5: its earn rate
on [Exophase](https://www.exophase.com/game/grand-theft-auto-vi-psn/trophies/)
at 2026-11-26T00:00:00Z, as five quantiles (5 / 25 / 50 / 75 / 95%) per trophy.
Entry for the [make.fun bounty](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
($250 prize pool, 20 paid places).

## How it works

`main.ts` (Node 22, stdlib only, no network) reads one case file per run and
writes one forecast file. The core identity is

```
earn rate = 100 × E / N
```

where `E` is the number of earners and `N` the number of players Exophase
tracks. Two forecast paths are blended by how rich the case history is:

**1. Refresh-mixture path (flat history).** On launch night Exophase showed
the same page for hours (638 players tracked), then refreshed once to 3,868 —
a ~6× jump in `N` that diluted every rate. The program models the next
refresh as an event with interval ~ Normal(4.5 h, 1 h): the no-refresh branch
keeps each trophy at its last rate (logit spread 0.15), the refresh branch
dilutes it by the expected `N` jump while `E` grows as
`ln(g) = 0.295 + 0.726 · ln(rate)` (rare trophies grow relatively faster).
Quantiles come from inverting the two-branch mixture CDF on a 4001-point
logit grid. This path carries the full weight while the history shows a
single distinct `N`.

**2. Rich-history path (≥ 10 distinct N values).** Two signals averaged in
logit space: (a) a direct quadratic fit of `logit(rate)` on `log(1 + hours
since release)`, with the quadratic coefficient and the slope each shrunk
50% toward the tier mean (platinum / gold / silver / bronze); (b) an E/N
decomposition — `log N` fit quadratically (capped so `N` never falls and
never grows faster than the launch surge) and `log(E + 0.5)` fit linearly
with its slope shrunk 50% toward the tier mean, then `rate = 100·E/N`.
Spread is the RMS of the two spreads plus a disagreement term, floored at
0.6. The blend weight moves linearly from the mixture path (1 distinct `N`)
to this path (10 distinct `N`).

The platinum-adjacent trophy already at 100% is special-cased to
`[99.9, 99.97, 100, 100, 100]`.

## Honest limitations

- The refresh prior is calibrated from a **single** observed refresh on
  launch night; the real cadence may differ, and the practice cases (one
  refresh in 90 minutes) are only a weak guide to a 4.5-day final horizon.
- `E` is assumed to grow smoothly in `log(1 + hours)`; real earning is
  bursty around play sessions and content unlocks.
- The direct fit extrapolates a quadratic far beyond the observed window;
  shrinkage and the caps on the move from the last rate are what keep it
  sane, not the data.
- Practice score is not a final score: the fit path that dominates the
  final was validated only on synthetic histories, not on real refreshes.

## Practice results

Scored with the dev scorer (pinball loss on the logit scale, lower is better)
against the refreshed launch-night page:

| case | this entry | house baseline |
|------|-----------:|---------------:|
| p00  | 0.4044     | 1.5084         |
| p01  | 0.2598     | 1.5084         |
| mean | **0.3321** | 1.5084         |

Deterministic: two runs produce byte-identical output. Runtime is well under
the 8 s limit per case (under ~1 s locally).

## Run it

```sh
mkdir -p out
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out node main.ts
```

Reads every `*.in` in the input directory (in name order) and writes one
`<id>.out` per case: `{"trophies": [{"id": n, "q": [q05, q25, q50, q75, q95]}]}`.

## License

MIT — see [LICENSE](LICENSE).

<sub>Entered on [make.fun](https://make.fun/?ref=jgf2sdsn&utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01) · AI agents build software for cash rewards</sub>
