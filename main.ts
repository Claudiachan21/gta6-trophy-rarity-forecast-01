// main.ts — perkiraan earn rate trofi GTA VI (PS5) di Exophase pada target case.
//
// Idenya: earn rate = 100 * E / N, dengan E = jumlah peraih dan N = pemain
// yang dilacak Exophase. Dua jalur perkiraan digabung berdasarkan seberapa
// kaya histori case-nya:
//
//   1. Jalur "campuran refresh" — dipakai penuh saat histori masih datar
//      (semua read menunjukkan N yang sama). Exophase me-refresh hitungannya
//      kira-kira tiap 4,5 jam; tiap refresh N melonjak beberapa kali lipat
//      sehingga semua rate terdilusi, sementara E ikut tumbuh mengikuti
//      pola ln(g) = 0,295 + 0,726 * ln(rate). Dua cabang (tanpa refresh vs
//      satu refresh sebelum target) dicampur lewat CDF di skala logit.
//
//   2. Jalur "histori kaya" — dipakai penuh saat histori sudah punya >= 10
//      nilai N berbeda. Dua sinyal dirata-rata di ruang logit: (a) fit
//      kuadratik langsung logit(rate) terhadap log(1+jam sejak rilis),
//      dengan koefisien kuadrat dan slope di-shrink 50% ke rata-rata tier;
//      (b) dekomposisi E/N: log N difit kuadratik, log(E+0,5) difit garis
//      dengan slope di-shrink 50% ke rata-rata tier, lalu rate = 100*E/N.
//
// Bobot campuran = 1 saat distinctN = 1, 0 saat distinctN >= 10, linear
// di antaranya. Trofi yang sudah 100% di-special-case.
//
// Node 22, stdlib saja, tanpa network, tanpa randomness: dua run pada case
// yang sama menghasilkan byte yang identik.

import { existsSync, mkdirSync, readdirSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";

// ---------------------------------------------------------------------------
// Konstanta
// ---------------------------------------------------------------------------

const RILIS_MS = Date.parse("2026-11-19T00:00:00Z"); // rilis GTA VI (UTC)
const LEVEL = [0.05, 0.25, 0.5, 0.75, 0.95];
const Z = [-1.645, -0.674, 0, 0.674, 1.645]; // kuantil normal standar
const GRID_N = 4001; // titik grid untuk inversi CDF campuran

// Prior refresh Exophase (dikalibrasi dari satu peristiwa: refresh malam
// peluncuran, N 638 -> 3868 dalam ~4,7 jam).
const REFRESH_RERATA_JAM = 4.5;
const REFRESH_SD_JAM = 1.0;
const N_TUMBUH_PER_JAM = 0.35; // log N per jam selama surge peluncuran
const LOMPAT_MAKS = 15; // cap pengali N saat refresh
const E_TUMBUH_A = 0.295; // ln(g) = A + B * ln(rate)
const E_TUMBUH_B = 0.726;
const E_TUMBUH_JAM_KALIBRASI = 4.7;
const SEBAR_TANPA_REFRESH = 0.15; // spread logit cabang tanpa refresh
const SEBAR_DENGAN_REFRESH = 0.5; // spread logit cabang refresh

// Jalur histori kaya.
const SEBAR_EN = 0.45; // spread tetap jalur E/N
const SEBAR_KAYA_LANTAI = 0.6;
const SEBAR_KAYA_ATAP = 1.4;

// ---------------------------------------------------------------------------
// Util matematika kecil
// ---------------------------------------------------------------------------

function jepit(x: number, lo: number, hi: number): number {
  return x < lo ? lo : x > hi ? hi : x;
}

function logit(persen: number): number {
  const p = jepit(persen, 0.01, 99.99) / 100;
  return Math.log(p / (1 - p));
}

function invLogit(l: number): number {
  return 100 / (1 + Math.exp(-l));
}

function bulat4(x: number): number {
  return Math.round(x * 10000) / 10000;
}

// Fungsi error (Abramowitz-Stegun 7.1.26, |galat| < 1,5e-7).
function erf(x: number): number {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const poly =
    (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
  const y = 1 - poly * Math.exp(-x * x);
  return x < 0 ? -y : y;
}

function phi(z: number): number {
  return 0.5 * (1 + erf(z / Math.SQRT2));
}

function rerata(xs: number[]): number {
  return xs.reduce((a, b) => a + b, 0) / xs.length;
}

function simpanganPopulasi(xs: number[], mu: number): number {
  if (xs.length < 2) return 0;
  const v = xs.reduce((a, x) => a + (x - mu) * (x - mu), 0) / xs.length;
  return Math.sqrt(v);
}

// Fit garis y = a + b x (kuadrat terkecil). Mengembalikan [a, b].
function pasangGaris(xs: number[], ys: number[]): [number, number] {
  const n = xs.length;
  const mx = rerata(xs);
  const my = rerata(ys);
  let sxx = 0;
  let sxy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (xs[i] - mx) * (xs[i] - mx);
    sxy += (xs[i] - mx) * (ys[i] - my);
  }
  const b = sxx > 1e-12 ? sxy / sxx : 0;
  return [my - b * mx, b];
}

// Fit kuadrat y = a + b x + c x^2 lewat persamaan normal + eliminasi Gauss
// dengan pivot parsial. Mengembalikan [a, b, c]; fallback ke garis datar
// bila sistemnya singular.
function pasangKuadrat(xs: number[], ys: number[]): [number, number, number] {
  const n = xs.length;
  const s = [0, 0, 0, 0, 0];
  const t = [0, 0, 0];
  for (let i = 0; i < n; i++) {
    const x = xs[i];
    const pw = [1, x, x * x, x * x * x, x * x * x * x];
    for (let k = 0; k < 5; k++) s[k] += pw[k];
    for (let k = 0; k < 3; k++) t[k] += pw[k] * ys[i];
  }
  const m = [
    [s[0], s[1], s[2], t[0]],
    [s[1], s[2], s[3], t[1]],
    [s[2], s[3], s[4], t[2]],
  ];
  for (let col = 0; col < 3; col++) {
    let piv = col;
    for (let r = col + 1; r < 3; r++) {
      if (Math.abs(m[r][col]) > Math.abs(m[piv][col])) piv = r;
    }
    const tmp = m[col];
    m[col] = m[piv];
    m[piv] = tmp;
    const d = m[col][col];
    if (Math.abs(d) < 1e-12) return [rerata(ys), 0, 0];
    for (let r = 0; r < 3; r++) {
      if (r === col) continue;
      const f = m[r][col] / d;
      for (let c = col; c < 4; c++) m[r][c] -= f * m[col][c];
    }
  }
  return [m[0][3] / m[0][0], m[1][3] / m[1][1], m[2][3] / m[2][2]];
}

// Jam sejak rilis (dijepit >= 0), lalu x = log(1 + jam) untuk fit.
function jamSejakRilis(iso: string): number {
  return Math.max(0, (Date.parse(iso) - RILIS_MS) / 3600000);
}

function kuantilDariMuSebar(mu: number, sebar: number): number[] {
  const qs = Z.map((z) => bulat4(jepit(invLogit(mu + z * sebar), 0, 100)));
  for (let k = 1; k < 5; k++) qs[k] = Math.max(qs[k], qs[k - 1]);
  return qs;
}

// ---------------------------------------------------------------------------
// Tipe data case (schema makefun/gta6-trophies-case/v1)
// ---------------------------------------------------------------------------

interface Trofi {
  id: number;
  name: string;
  tier: string;
}

interface Read {
  at: string;
  players_tracked: number | null;
  rates: (number | null)[];
}

interface Kasus {
  id: string;
  at: string;
  horizon_minutes: number;
  target: string;
  trophies: Trofi[];
  history: Read[];
}

// Baris spesial untuk trofi yang sudah 100%: praktis tidak bisa turun.
const BARIS_SERATUS = [99.9, 99.97, 100, 100, 100];

// ---------------------------------------------------------------------------
// Jalur 1: campuran refresh (histori datar)
// ---------------------------------------------------------------------------

function jalurCampuran(k: Kasus): number[][] {
  const reads = k.history;
  const akhir = reads[reads.length - 1];
  const jamHorizon = k.horizon_minutes / 60;
  // S = berapa lama halaman sudah tidak berubah saat checkpoint.
  const basiJam = reads.length > 0
    ? Math.max(0, (Date.parse(k.at) - Date.parse(reads[0].at)) / 3600000)
    : 0;
  // Peluang setidaknya satu refresh sebelum target: interval refresh
  // dimodelkan Normal(4,5 jam, 1 jam).
  const pRefresh = phi((basiJam + jamHorizon - REFRESH_RERATA_JAM) / REFRESH_SD_JAM);
  // Lompatan N bila refresh: tumbuh log selama (basi + horizon) jam.
  const lompat = Math.min(Math.exp(N_TUMBUH_PER_JAM * (basiJam + jamHorizon)), LOMPAT_MAKS);
  const totalJam = basiJam + jamHorizon;
  const skalaTumbuh = totalJam / E_TUMBUH_JAM_KALIBRASI;

  // Grid logit untuk inversi CDF campuran (dibangun sekali).
  const lo = logit(0.01);
  const hi = logit(99.99);
  const grid: number[] = [];
  for (let i = 0; i < GRID_N; i++) grid.push(lo + ((hi - lo) * i) / (GRID_N - 1));

  const hasil: number[][] = [];
  for (let i = 0; i < k.trophies.length; i++) {
    let r0 = akhir.rates[i];
    if (r0 === null || r0 === undefined) r0 = 0.01;
    if (r0 >= 99.99) {
      hasil.push([...BARIS_SERATUS]);
      continue;
    }
    // Pertumbuhan E selama interval refresh, bergantung pada rate awal:
    // trofi langka tumbuh lebih cepat (relatif) daripada trofi umum.
    let g = Math.exp((E_TUMBUH_A + E_TUMBUH_B * Math.log(Math.max(r0, 0.01))) * skalaTumbuh);
    g = jepit(g, 0.8, 8);
    const rCair = Math.max(0.01, (r0 * g) / lompat);
    const muDiam = logit(r0);
    const muCair = logit(rCair);
    // CDF campuran di grid, lalu balikkan pada tiap level.
    const cdf = grid.map(
      (x) =>
        (1 - pRefresh) * phi((x - muDiam) / SEBAR_TANPA_REFRESH) +
        pRefresh * phi((x - muCair) / SEBAR_DENGAN_REFRESH),
    );
    const qs: number[] = [];
    for (const tau of LEVEL) {
      let idx = 0;
      while (idx < GRID_N && cdf[idx] < tau) idx++;
      if (idx >= GRID_N) idx = GRID_N - 1;
      qs.push(bulat4(jepit(invLogit(grid[idx]), 0, 100)));
    }
    for (let j = 1; j < 5; j++) qs[j] = Math.max(qs[j], qs[j - 1]);
    hasil.push(qs);
  }
  return hasil;
}

// ---------------------------------------------------------------------------
// Jalur 2: histori kaya (fit langsung + dekomposisi E/N)
// ---------------------------------------------------------------------------

interface Sinyal {
  mu: number;
  sebar: number;
}

// Slope rata-rata per tier dari fit garis tiap trofi (dipakai untuk shrink).
function slopeRataTier(
  k: Kasus,
  xsPerTrofi: number[][],
  ysPerTrofi: number[][],
): Record<string, number> {
  const out: Record<string, number> = {};
  for (const tier of ["platinum", "gold", "silver", "bronze"]) {
    const slopes: number[] = [];
    for (let i = 0; i < k.trophies.length; i++) {
      if (k.trophies[i].tier !== tier) continue;
      if (xsPerTrofi[i].length >= 3) slopes.push(pasangGaris(xsPerTrofi[i], ysPerTrofi[i])[1]);
    }
    out[tier] = slopes.length > 0 ? rerata(slopes) : 0;
  }
  return out;
}

function jalurHistori(k: Kasus): number[][] {
  const reads = k.history;
  const akhir = reads[reads.length - 1];
  const hariHorizon = k.horizon_minutes / 1440;
  const xTarget = Math.log(1 + jamSejakRilis(k.target));

  // Absis tiap read: log(1 + jam sejak rilis).
  const xs = reads.map((r) => Math.log(1 + jamSejakRilis(r.at)));

  // --- Sinyal A: fit kuadratik langsung pada logit(rate). ---
  const xsA: number[][] = [];
  const ysA: number[][] = [];
  for (let i = 0; i < k.trophies.length; i++) {
    const px: number[] = [];
    const py: number[] = [];
    for (let j = 0; j < reads.length; j++) {
      const r = reads[j].rates[i];
      if (r !== null && r !== undefined) {
        px.push(xs[j]);
        py.push(logit(r));
      }
    }
    xsA.push(px);
    ysA.push(py);
  }
  const slopeTierA = slopeRataTier(k, xsA, ysA);

  // --- Sinyal B: dekomposisi E/N. ---
  // Log N difit kuadratik terhadap x.
  const idxN: number[] = [];
  const logN: number[] = [];
  for (let j = 0; j < reads.length; j++) {
    const n = reads[j].players_tracked;
    if (n !== null && n !== undefined && n > 0) {
      idxN.push(xs[j]);
      logN.push(Math.log(n));
    }
  }
  const nAkhir = akhir.players_tracked && akhir.players_tracked > 0 ? akhir.players_tracked : 1;
  let logNT: number;
  if (idxN.length >= 3) {
    const [a, b, c] = pasangKuadrat(idxN, logN);
    logNT = a + b * xTarget + c * xTarget * xTarget;
  } else {
    logNT = Math.log(nAkhir);
  }
  // N tak boleh turun, dan tak boleh tumbuh lebih cepat dari surge awal.
  logNT = jepit(
    logNT,
    Math.log(nAkhir),
    Math.log(nAkhir) + 0.5 * Math.sqrt(hariHorizon) + 0.35 * hariHorizon,
  );

  const xsB: number[][] = [];
  const ysB: number[][] = [];
  for (let i = 0; i < k.trophies.length; i++) {
    const px: number[] = [];
    const py: number[] = [];
    for (let j = 0; j < reads.length; j++) {
      const r = reads[j].rates[i];
      const n = reads[j].players_tracked;
      if (r !== null && r !== undefined && n !== null && n !== undefined && n > 0) {
        px.push(xs[j]);
        py.push(Math.log((r * n) / 100 + 0.5)); // log(E + 0,5)
      }
    }
    xsB.push(px);
    ysB.push(py);
  }
  const slopeTierB = slopeRataTier(k, xsB, ysB);

  // --- Gabungkan kedua sinyal per trofi. ---
  const hasil: number[][] = [];
  for (let i = 0; i < k.trophies.length; i++) {
    let r0 = akhir.rates[i];
    if (r0 === null || r0 === undefined) r0 = 0.01;
    if (r0 >= 99.99) {
      hasil.push([...BARIS_SERATUS]);
      continue;
    }
    const tier = k.trophies[i].tier;

    // Sinyal A.
    let sinyalA: Sinyal;
    if (xsA[i].length >= 3) {
      let [a, b, c] = pasangKuadrat(xsA[i], ysA[i]);
      c *= 0.5; // shrink kelengkungan 50%
      b = 0.5 * b + 0.5 * slopeTierA[tier]; // shrink slope 50% ke mean tier
      // Jangkar ulang intercept pada rata-rata agar shrink tak geser level.
      const mx = rerata(xsA[i]);
      const mx2 = rerata(xsA[i].map((x) => x * x));
      a = rerata(ysA[i]) - b * mx - c * mx2;
      let mu = a + b * xTarget + c * xTarget * xTarget;
      mu = jepit(mu, logit(r0) - 3.5, logit(r0) + 1.0);
      const resid = xsA[i].map((x, j) => ysA[i][j] - (a + b * x + c * x * x));
      const rsd = simpanganPopulasi(resid, rerata(resid));
      const sebar = jepit(
        Math.sqrt(
          Math.max(0.15, rsd) ** 2 +
            (0.2 * Math.sqrt(hariHorizon)) ** 2 +
            (0.25 * Math.sqrt(hariHorizon)) ** 2,
        ),
        0.3,
        1.0,
      );
      sinyalA = { mu, sebar };
    } else {
      sinyalA = { mu: logit(r0), sebar: 0.5 };
    }

    // Sinyal B.
    let sinyalB: Sinyal;
    if (xsB[i].length >= 3) {
      let [a, b] = pasangGaris(xsB[i], ysB[i]);
      b = 0.5 * b + 0.5 * slopeTierB[tier];
      a = rerata(ysB[i]) - b * rerata(xsB[i]);
      const logET = Math.max(a + b * xTarget, Math.log(0.5));
      const rateT = jepit(100 * Math.exp(logET - logNT), 0.01, 99.99);
      sinyalB = { mu: logit(rateT), sebar: SEBAR_EN };
    } else {
      sinyalB = { mu: logit(r0), sebar: 0.5 };
    }

    // Rata-rata median di ruang logit; spread = RMS kedua spread ditambah
    // suku ketidaksepakatan, dengan lantai 0,6.
    const mu = 0.5 * (sinyalA.mu + sinyalB.mu);
    const selisih = sinyalA.mu - sinyalB.mu;
    const sebar = jepit(
      Math.sqrt(0.5 * sinyalA.sebar ** 2 + 0.5 * sinyalB.sebar ** 2 + 0.25 * selisih * selisih),
      SEBAR_KAYA_LANTAI,
      SEBAR_KAYA_ATAP,
    );
    hasil.push(kuantilDariMuSebar(mu, sebar));
  }
  return hasil;
}

// ---------------------------------------------------------------------------
// Perkiraan satu case + main
// ---------------------------------------------------------------------------

function ramalKasus(k: Kasus): { id: number; q: number[] }[] {
  // Berapa nilai N berbeda yang pernah terlihat: ukuran kekayaan histori.
  const nilaiN = new Set<number>();
  for (const r of k.history) {
    if (r.players_tracked !== null && r.players_tracked !== undefined) {
      nilaiN.add(r.players_tracked);
    }
  }
  const distinctN = nilaiN.size;
  // Bobot jalur campuran: 1 saat histori datar, 0 saat histori kaya.
  const bobotCampuran = jepit((10 - distinctN) / 9, 0, 1);

  const campuran = jalurCampuran(k);
  const histori = jalurHistori(k);

  const baris: { id: number; q: number[] }[] = [];
  for (let i = 0; i < k.trophies.length; i++) {
    let q: number[];
    if (bobotCampuran >= 1) {
      q = campuran[i];
    } else if (bobotCampuran <= 0) {
      q = histori[i];
    } else {
      // Blend tiap kuantil di ruang logit (tetap non-decreasing).
      q = campuran[i].map((qc, j) =>
        bulat4(jepit(invLogit(bobotCampuran * logit(qc) + (1 - bobotCampuran) * logit(histori[i][j])), 0, 100)),
      );
      for (let j = 1; j < 5; j++) q[j] = Math.max(q[j], q[j - 1]);
    }
    baris.push({ id: k.trophies[i].id, q });
  }
  return baris;
}

function main(): void {
  const dirMasuk = process.env.MAKEFUN_INPUT_DIR || "/input/";
  const dirKeluar = process.env.MAKEFUN_OUTPUT_DIR || "/output/";
  if (!existsSync(dirKeluar)) mkdirSync(dirKeluar, { recursive: true });
  const files = readdirSync(dirMasuk)
    .filter((f) => f.endsWith(".in"))
    .sort();
  for (const f of files) {
    const kasus: Kasus = JSON.parse(readFileSync(join(dirMasuk, f), "utf-8"));
    const trofi = ramalKasus(kasus);
    const namaKeluar = f.replace(/\.in$/, ".out");
    writeFileSync(join(dirKeluar, namaKeluar), JSON.stringify({ trophies: trofi }) + "\n", "utf-8");
  }
}

main();
