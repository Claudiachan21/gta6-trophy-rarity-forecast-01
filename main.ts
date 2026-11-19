// GTA VI trophy rarity forecast: a dilution-aware forecaster.
//
// The thing being forecast is a ratio: each trophy's earn rate on Exophase is
//   rate = (players who earned the trophy) / (players Exophase tracks).
// Both the numerator and the denominator grow after launch, but they grow at
// very different speeds, and that is the whole game. On launch night the
// tracked-player count sat at 638 for hours while the page went stale, then
// one refresh took it to 3,868 and every rate but the one at 100% fell: the
// denominator jumped 6x while the numerators barely moved. New players join
// far faster than they earn trophies, so rates dilute down, hardest for the
// trophies nobody earns in their first hours and least for the ones every new
// player earns at once (the 100% trophy never moves at all).
//
// The program therefore never forecasts a rate directly. It forecasts the two
// counts separately and divides:
//   - N(t), the tracked players, from a fit of log N against log time;
//   - E_i(t), the earners of trophy i, from a fit of log E_i against log time
//     (shrunk toward its tier's slope, because a single trophy's history is
//     noisy when its count is small);
//   - rate_i = 100 * E_i / N.
// When the history is flat (every read identical, the page stale), no trend
// can be fit, so the program falls back to a refresh prior: Exophase
// refreshes roughly every few hours, and the longer the page has sat
// unchanged the likelier a refresh lands inside the horizon. The forecast
// then becomes a mixture of "no refresh, rate stays" and "refresh, rate
// dilutes by the expected jump", with the mixture weight from the staleness.
//
// The quantiles are honest: the spread combines the fit residuals, the
// uncertainty of extrapolating N, and a random-walk term that grows with the
// horizon. The pinball loss is a proper scoring rule, so the program writes
// the quantiles it believes, never a sharpened or hedged version of them.
//
// Everything here is deterministic: no clock, no randomness, no network. Two
// runs on the same case write the same bytes.

import { mkdirSync, readdirSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";

const INPUT_DIR = process.env.MAKEFUN_INPUT_DIR ?? "/input";
const OUTPUT_DIR = process.env.MAKEFUN_OUTPUT_DIR ?? "/output";

// The game released at the start of 2026-11-19 UTC; log time is measured
// from there, so the launch-night reads sit near zero.
const RELEASE_MS = Date.parse("2026-11-19T00:00:00Z");
const HOUR_MS = 3_600_000;

// The five quantile levels the contest scores, in the order each row lists.
const LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95];
// The standard normal's quantile at each level.
const Z = [-1.645, -0.674, 0, 0.674, 1.645];

// --- small numeric helpers -------------------------------------------------

function clampRate(pct: number): number {
  return Math.min(99.99, Math.max(0.01, pct));
}

function logit(pct: number): number {
  const p = clampRate(pct) / 100;
  return Math.log(p / (1 - p));
}

function invLogit(l: number): number {
  return 100 / (1 + Math.exp(-l));
}

function round4(x: number): number {
  return Math.round(x * 10_000) / 10_000;
}

function erfApprox(x: number): number {
  // Abramowitz and Stegun 7.1.26: |error| < 1.5e-7, plenty for a grid CDF.
  const s = x < 0 ? -1 : 1;
  const ax = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * ax);
  const poly =
    (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
  return s * (1 - poly * Math.exp(-ax * ax));
}

function normCdf(z: number): number {
  return 0.5 * (1 + erfApprox(z / Math.SQRT2));
}

function mean(xs: number[]): number {
  return xs.reduce((a, b) => a + b, 0) / xs.length;
}

function sd(xs: number[]): number {
  if (xs.length < 2) return 0;
  const m = mean(xs);
  return Math.sqrt(xs.reduce((a, b) => a + (b - m) * (b - m), 0) / (xs.length - 1));
}

// Least squares for y = a + b*x. Returns [a, b].
function fitLine(xs: number[], ys: number[]): [number, number] {
  const n = xs.length;
  const mx = mean(xs);
  const my = mean(ys);
  let sxx = 0;
  let sxy = 0;
  for (let k = 0; k < n; k += 1) {
    sxx += (xs[k] - mx) * (xs[k] - mx);
    sxy += (xs[k] - mx) * (ys[k] - my);
  }
  const b = sxx > 1e-12 ? sxy / sxx : 0;
  return [my - b * mx, b];
}

// Least squares for y = a + b*x + c*x^2, solved by plain Gaussian
// elimination on the 3x3 normal equations. Returns [a, b, c].
function fitQuad(xs: number[], ys: number[]): [number, number, number] {
  const n = xs.length;
  let s0 = n;
  let s1 = 0;
  let s2 = 0;
  let s3 = 0;
  let s4 = 0;
  let t0 = 0;
  let t1 = 0;
  let t2 = 0;
  for (let k = 0; k < n; k += 1) {
    const x = xs[k];
    const y = ys[k];
    const x2 = x * x;
    s1 += x;
    s2 += x2;
    s3 += x2 * x;
    s4 += x2 * x2;
    t0 += y;
    t1 += x * y;
    t2 += x2 * y;
  }
  // The system [[s0,s1,s2],[s1,s2,s3],[s2,s3,s4]] * [a,b,c] = [t0,t1,t2].
  const m = [
    [s0, s1, s2, t0],
    [s1, s2, s3, t1],
    [s2, s3, s4, t2],
  ];
  for (let col = 0; col < 3; col += 1) {
    let piv = col;
    for (let r = col + 1; r < 3; r += 1) {
      if (Math.abs(m[r][col]) > Math.abs(m[piv][col])) piv = r;
    }
    const tmp = m[col];
    m[col] = m[piv];
    m[piv] = tmp;
    const d = m[col][col];
    if (Math.abs(d) < 1e-12) return [t0 / Math.max(n, 1), 0, 0];
    for (let r = 0; r < 3; r += 1) {
      if (r === col) continue;
      const f = m[r][col] / d;
      for (let c = col; c < 4; c += 1) m[r][c] -= f * m[col][c];
    }
  }
  return [m[0][3] / m[0][0], m[1][3] / m[1][1], m[2][3] / m[2][2]];
}

// --- case shape ------------------------------------------------------------

type Trophy = { id: number; name: string; tier: string };
type Read = { at: string; players_tracked: number; rates: Array<number | null> };
type CaseFile = {
  id: string;
  at: string;
  horizon_minutes: number;
  target: string;
  levels: number[];
  trophies: Trophy[];
  history: Read[];
};

// --- the refresh prior (used when the history is flat) ----------------------
// These numbers are the model, fitted once from the launch-night refresh the
// rules describe (the page sat unchanged at 638 players, then one refresh
// took it to 3,868 and every rate but the 100% one fell) and from the
// earn-count growth that refresh revealed. They are priors, not the answers:
// the program applies them to whatever case it is given, and the fit takes
// over whenever the history has anything to fit.
const REFRESH_MEAN_H = 4.5; // Exophase's refresh cadence, in hours
const REFRESH_SD_H = 1.0;
const N_GROWTH_PER_H = 0.35; // log growth of tracked players in the surge
const E_GROWTH_A = 0.295; // ln(E growth) = A + B * ln(rate): new players earn
const E_GROWTH_B = 0.726; // the early, high-rate trophies fastest
const E_GROWTH_SPAN_H = 4.7; // the span that E-growth fit was calibrated over

// The quantiles of a two-branch mixture (no refresh vs refresh), computed by
// walking a fine grid over the logit scale and inverting the mixture CDF.
// Deterministic: the grid is fixed, so two runs agree byte for byte.
function mixtureQuantiles(
  rateFlat: number,
  rateDiluted: number,
  pRefresh: number,
  spreadFlat: number,
  spreadDiluted: number,
): number[] {
  const lo = logit(0.01);
  const hi = logit(99.99);
  const steps = 4000;
  const grid: number[] = [];
  for (let k = 0; k <= steps; k += 1) grid.push(lo + ((hi - lo) * k) / steps);
  const muFlat = logit(rateFlat);
  const muDiluted = logit(rateDiluted);
  const cdf: number[] = [];
  for (let k = 0; k <= steps; k += 1) {
    const x = grid[k];
    cdf.push(
      (1 - pRefresh) * normCdf((x - muFlat) / spreadFlat) +
        pRefresh * normCdf((x - muDiluted) / spreadDiluted),
    );
  }
  const qs: number[] = [];
  for (const tau of LEVELS) {
    let j = 0;
    while (j <= steps && cdf[j] < tau) j += 1;
    if (j > steps) j = steps;
    qs.push(round4(invLogit(grid[j])));
  }
  for (let k = 1; k < qs.length; k += 1) qs[k] = Math.max(qs[k], qs[k - 1]);
  return qs;
}

// --- the direct rate fit: quadratic logit(rate) on log(1+hours) ---------
// The second rich-history signal. It is exact when dilution is a power law
// in time; the E/N split is better when the player base saturates. The two
// disagree on different worlds, so the rich forecast averages them.
function directFit(
  reads: Read[],
  xs: number[],
  xTarget: number,
  trophyIndex: number,
  tier: string,
  tierSlope: Map<string, number>,
  horizonDays: number,
  lastRate: number,
): { median: number; spread: number } {
  const pts: Array<[number, number]> = [];
  reads.forEach((r, j) => {
    const rate = r.rates[trophyIndex];
    if (rate === null || rate === undefined) return;
    pts.push([xs[j], logit(rate)]);
  });
  if (pts.length < 3) return { median: lastRate, spread: 0.5 };
  const q = fitQuad(
    pts.map((p) => p[0]),
    pts.map((p) => p[1]),
  );
  // Shrink the curvature halfway toward flat and the slope halfway toward
  // the tier mean; re-anchor the intercept at the data mean so shrinkage
  // does not shift the level.
  const c = 0.5 * q[2];
  const meanX = mean(pts.map((p) => p[0]));
  const meanX2 = mean(pts.map((p) => p[0] * p[0]));
  const meanY = mean(pts.map((p) => p[1]));
  const b = 0.5 * q[1] + 0.5 * (tierSlope.get(tier) ?? 0);
  const a = meanY - b * meanX - c * meanX2;
  let mu = a + b * xTarget + c * xTarget * xTarget;
  // Never move more than 3.5 logit down / 1.0 up from the last read rate.
  mu = Math.max(logit(lastRate) - 3.5, Math.min(logit(lastRate) + 1.0, mu));
  const residSd = Math.sqrt(
    Math.max(
      0.0225,
      mean(pts.map(([x, y]) => (y - (a + b * x + c * x * x)) ** 2)),
    ),
  );
  const spread = Math.max(
    0.3,
    Math.min(
      1.0,
      Math.sqrt(
        residSd * residSd +
          (0.2 * Math.sqrt(horizonDays)) ** 2 +
          (0.25 * Math.sqrt(horizonDays)) ** 2,
      ),
    ),
  );
  return { median: clampRate(invLogit(mu)), spread };
}

// --- the forecast ----------------------------------------------------------

function forecastCase(c: CaseFile): Array<{ id: number; q: number[] }> {
  const reads = c.history;
  if (reads.length === 0) throw new Error(`${c.id} has no history`);
  const last = reads[reads.length - 1];
  const horizonDays = c.horizon_minutes / 1440;
  const horizonH = c.horizon_minutes / 60;

  // Hours since release for each read and for the target, clamped at zero
  // (reads from just before midnight sit at the release itself).
  const hoursOf = (iso: string): number =>
    Math.max(0, (Date.parse(iso) - RELEASE_MS) / HOUR_MS);
  const readHours = reads.map((r) => hoursOf(r.at));
  const targetHours = hoursOf(c.target);
  const xOf = (h: number): number => Math.log(1 + h);
  const xs = readHours.map(xOf);
  const xTarget = xOf(targetHours);

  // How many distinct player counts the history shows. One means the page
  // sat unchanged the whole time (the practice cases); many means refreshes
  // happened and a trend can be fit.
  const distinctN = new Set(reads.map((r) => r.players_tracked)).size;

  // --- the refresh prior, always computed (it is the whole forecast when
  // the history is flat, and a blend weight when it is not) ---
  const stalenessH = Math.max(
    0,
    (Date.parse(c.at) - Date.parse(reads[0].at)) / HOUR_MS,
  );
  const pRefresh = normCdf((stalenessH + horizonH - REFRESH_MEAN_H) / REFRESH_SD_H);
  const jump = Math.min(15, Math.exp(N_GROWTH_PER_H * (stalenessH + horizonH)));
  const eSpan = Math.max(0.5, (stalenessH + horizonH) / E_GROWTH_SPAN_H);

  // --- the N fit (log players against log time, quadratic) ---
  const logN = reads.map((r) => Math.log(Math.max(1, r.players_tracked)));
  const [aN, bN, cN] = distinctN > 1 ? fitQuad(xs, logN) : [logN[logN.length - 1], 0, 0];
  let logNTarget = aN + bN * xTarget + cN * xTarget * xTarget;
  // Players only join: the forecast never goes below the last count, and it
  // never grows faster than the launch surge sustained.
  logNTarget = Math.max(logN[logN.length - 1], logNTarget);
  logNTarget = Math.min(
    logNTarget,
    logN[logN.length - 1] + 0.5 * Math.sqrt(horizonDays) + 0.35 * horizonDays,
  );
  const nTarget = Math.exp(logNTarget);

  // --- the per-tier earner slopes, for shrinkage ---
  const tiers = ["platinum", "gold", "silver", "bronze"];
  const tierSlope = new Map<string, number>();
  const tierRateSlope = new Map<string, number>();
  for (const tier of tiers) {
    const slopes: number[] = [];
    const rateSlopes: number[] = [];
    c.trophies.forEach((t, i) => {
      if (t.tier !== tier) return;
      const pts: Array<[number, number]> = [];
      const ratePts: Array<[number, number]> = [];
      reads.forEach((r, j) => {
        const rate = r.rates[i];
        if (rate === null || rate === undefined) return;
        const e = Math.max(0, (rate * r.players_tracked) / 100);
        pts.push([xs[j], Math.log(e + 0.5)]);
        ratePts.push([xs[j], logit(rate)]);
      });
      if (pts.length >= 3) {
        const [, b] = fitLine(
          pts.map((p) => p[0]),
          pts.map((p) => p[1]),
        );
        slopes.push(b);
        const [, br] = fitLine(
          ratePts.map((p) => p[0]),
          ratePts.map((p) => p[1]),
        );
        rateSlopes.push(br);
      }
    });
    tierSlope.set(tier, slopes.length > 0 ? mean(slopes) : 0);
    tierRateSlope.set(tier, rateSlopes.length > 0 ? mean(rateSlopes) : 0);
  }

  // --- per trophy ---
  const rows: Array<{ id: number; q: number[] }> = [];
  c.trophies.forEach((t, i) => {
    const r0 = last.rates[i] === null || last.rates[i] === undefined ? 0.01 : last.rates[i] as number;

    // The trophy everyone has: it stays at 100% whatever the denominator
    // does, because every new player earns it too.
    if (r0 >= 99.99) {
      rows.push({ id: t.id, q: [99.9, 99.97, 100, 100, 100] });
      return;
    }

    // --- the fit forecast: earners over players ---
    const pts: Array<[number, number]> = [];
    reads.forEach((r, j) => {
      const rate = r.rates[i];
      if (rate === null || rate === undefined) return;
      const e = Math.max(0, (rate * r.players_tracked) / 100);
      pts.push([xs[j], Math.log(e + 0.5)]);
    });
    let fitMedian = r0;
    let fitSpread = 0.5;
    if (pts.length >= 3 && distinctN > 1) {
      const [aE, bE] = fitLine(
        pts.map((p) => p[0]),
        pts.map((p) => p[1]),
      );
      // Shrink the trophy's slope halfway toward its tier's: one trophy's
      // history is noisy when its count is small, and trophies of a tier
      // are earned at similar speeds.
      const bTier = tierSlope.get(t.tier) ?? 0;
      const bShrunk = 0.5 * bE + 0.5 * bTier;
      const aShrunk = mean(pts.map((p) => p[1])) - bShrunk * mean(pts.map((p) => p[0]));
      let logETarget = aShrunk + bShrunk * xTarget;
      // Earners never fall, and never exceed the players there are.
      const eLast = Math.max(0, (r0 * last.players_tracked) / 100);
      logETarget = Math.max(Math.log(eLast + 0.5), logETarget);
      logETarget = Math.min(logETarget, Math.log(nTarget + 0.5));
      const eTarget = Math.exp(logETarget) - 0.5;
      fitMedian = clampRate((100 * Math.max(0, eTarget)) / nTarget);
      const resid = pts.map(([x, y]) => y - (aShrunk + bShrunk * x));
      const residSd = sd(resid);
      fitSpread = Math.min(
        1.0,
        Math.sqrt(
          Math.max(0.15, residSd) ** 2 +
            (0.2 * Math.sqrt(horizonDays)) ** 2 +
            (0.25 * Math.sqrt(horizonDays)) ** 2,
        ),
      );
      fitSpread = Math.max(0.3, fitSpread);
    }
    // --- the rich-history forecast: average of the two fit signals ---
    // The direct rate fit and the E/N split are exact on different worlds
    // (power-law dilution vs a saturating player base), so average their
    // medians in logit space. The spread is the root-mean-square of the two
    // spreads plus a disagreement term, floored at 0.6: the two methods
    // disagreeing is itself evidence the future is uncertain.
    const direct = directFit(reads, xs, xTarget, i, t.tier, tierRateSlope, horizonDays, r0);
    const muDirect = logit(direct.median);
    const muEn = logit(fitMedian);
    const richMu = 0.5 * muDirect + 0.5 * muEn;
    const richSpread = Math.max(
      0.6,
      Math.min(
        1.4,
        Math.sqrt(
          0.5 * direct.spread * direct.spread +
            0.5 * fitSpread * fitSpread +
            0.25 * (muDirect - muEn) * (muDirect - muEn),
        ),
      ),
    );
    const fitQs = Z.map((z) => round4(invLogit(richMu + z * richSpread)));

    // --- the refresh-mixture forecast (the flat-history answer) ---
    const gRaw = Math.exp((E_GROWTH_A + E_GROWTH_B * Math.log(Math.max(r0, 0.01))) * eSpan);
    const g = Math.min(8, Math.max(0.8, gRaw));
    const rDiluted = clampRate((r0 * g) / jump);
    const mixQs = mixtureQuantiles(r0, rDiluted, pRefresh, 0.15, 0.5);

    // --- blend: flat history is all mixture, rich history is all fit ---
    let wMix: number;
    if (distinctN <= 1) wMix = 1;
    else if (distinctN >= 10) wMix = 0;
    else wMix = (10 - distinctN) / 9;
    const qs = LEVELS.map((_, k) => {
      const blended = invLogit(
        wMix * logit(mixQs[k]) + (1 - wMix) * logit(fitQs[k]),
      );
      return round4(clampRate(blended));
    });
    for (let k = 1; k < qs.length; k += 1) qs[k] = Math.max(qs[k], qs[k - 1]);
    rows.push({ id: t.id, q: qs });
  });
  return rows;
}

// --- driver ----------------------------------------------------------------

mkdirSync(OUTPUT_DIR, { recursive: true });
const inputs = readdirSync(INPUT_DIR)
  .filter((f) => f.endsWith(".in"))
  .sort();
if (inputs.length === 0) throw new Error(`no .in files in ${INPUT_DIR}`);
for (const file of inputs) {
  const id = file.slice(0, -".in".length);
  const c = JSON.parse(readFileSync(join(INPUT_DIR, file), "utf8")) as CaseFile;
  if (
    c.levels.length !== LEVELS.length ||
    c.levels.some((t, j) => t !== LEVELS[j])
  ) {
    throw new Error(`${file} asks for levels ${JSON.stringify(c.levels)}`);
  }
  const trophies = forecastCase(c);
  if (trophies.length !== c.trophies.length) {
    throw new Error(`${file}: forecast ${trophies.length} rows, case has ${c.trophies.length}`);
  }
  const out = `${JSON.stringify({ trophies })}\n`;
  writeFileSync(join(OUTPUT_DIR, `${id}.out`), out);
}
