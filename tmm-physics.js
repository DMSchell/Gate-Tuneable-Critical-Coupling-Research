// MARK: Complex Functions
// -----------------------------------------------------------------------------------------------

function cAdd(a, b) { return { re: a.re + b.re, im: a.im + b.im }; }
function cSub(a, b) { return { re: a.re - b.re, im: a.im - b.im }; }
function cMul(a, b) { return { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }; }

function cDiv(a, b) {
  const denom = b.re * b.re + b.im * b.im;
  return { re: (a.re * b.re + a.im * b.im) / denom, im: (a.im * b.re - a.re * b.im) / denom };
}

function cSqrt(a) {
  const r = Math.hypot(a.re, a.im);
  const re = Math.sqrt((r + a.re) / 2);
  let im = Math.sqrt((r - a.re) / 2);
  if (a.im < 0) im = -im;
  return { re, im };
}

function cExp(a) {
  const mag = Math.exp(a.re);
  return { re: mag * Math.cos(a.im), im: mag * Math.sin(a.im) };
}

function cScaleR(a, s) { return { re: a.re * s, im: a.im * s }; }
function cAbs2(a) { return a.re * a.re + a.im * a.im; }

// MARK: Constants + Parameters
// -----------------------------------------------------------------------------------------------

const E_CHARGE = 1.602176634e-19;
const HBAR     = 1.054571817e-34;
const EPS0     = 8.8541878128e-12;
const KB_EV    = 8.617333262e-5;
const ETA0     = 376.730313668;
const EV_NM    = 1239.841984;

const E2_OVER_HBAR = (E_CHARGE * E_CHARGE) / HBAR;
const SIGMA0       = E2_OVER_HBAR / 4;

const GRAPHENE_DEFAULTS = {
  VDirac:   0,
  epsOx:    3.9,
  tox_nm:   285,
  T_K:      300,
  vF:       1.0e6,
  nimp_cm2: 1e11,
  epsG:     2.45,
};

// MARK: Gate Electrostatics
// -----------------------------------------------------------------------------------------------
function grapheneGate(Vg, gp) {
  const Cg = (EPS0 * gp.epsOx) / (gp.tox_nm * 1e-9);
  const n_m2 = (Cg * (Vg - gp.VDirac)) / E_CHARGE;
  const mu_J = HBAR * gp.vF * Math.sqrt(Math.PI * Math.abs(n_m2)) * Math.sign(n_m2);
  return { n_m2, mu_eV: mu_J / E_CHARGE };
}

function grapheneScattering(mu_eV, T_eV, gp) {
  const echar_J = Math.max(Math.abs(mu_eV), T_eV, 1e-3) * E_CHARGE;
  const nimp_m2 = gp.nimp_cm2 * 1e4;
  const invTau = Math.pow(E_CHARGE, 4) * nimp_m2 / (8 * EPS0 * EPS0 * HBAR * gp.epsG * gp.epsG * echar_J); // 1/s
  return { tau_s: 1 / invTau, Gamma_eV: (HBAR * invTau) / E_CHARGE };
}

// MARK: Kubo Conductivity
// -----------------------------------------------------------------------------------------------
function ln2cosh(x) {
  const ax = Math.abs(x);
  return ax + Math.log1p(Math.exp(-2 * ax));
}

function grapheneG(x, mu, T) {
  const u = x / T, m = mu / T, M = Math.max(u, m);
  const eu = Math.exp(u - M), eun = Math.exp(-u - M);
  const em = Math.exp(m - M), emn = Math.exp(-m - M);
  return (eu - eun) / (em + emn + eu + eun);
}

function grapheneGprime(x, mu, T) {
  const u = x / T, m = mu / T, M = Math.max(u, m);
  const a = Math.exp(u - M) + Math.exp(-u - M);
  const b = Math.exp(m - M) + Math.exp(-m - M);
  const num = 0.25 * a * b + Math.exp(-2 * M);
  const den = T * Math.pow((a + b) / 2, 2);
  return num / den;
}

function simpson(f, a, b, n) {
  if (b - a <= 0) return 0;
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let k = 1; k < n; k++) s += f(a + k * h) * (k % 2 ? 4 : 2);
  return (s * h) / 3;
}

function grapheneSigma(E_eV, mu_eV, T_eV, Gamma_eV) {
  const w  = Math.max(E_eV, 1e-6);
  const T  = Math.max(T_eV, 1e-5);
  const mu = Math.abs(mu_eV);
  const G  = Math.max(Gamma_eV, 1e-9);

  const pref = (2 * E2_OVER_HBAR / Math.PI) * T * ln2cosh(mu / (2 * T));
  const d2 = w * w + G * G;
  let sRe = (pref * G) / d2;
  let sIm = (pref * w) / d2;

  const G0  = grapheneG(w / 2, mu, T);
  const Gp0 = grapheneGprime(w / 2, mu, T);
  const Lam = Math.max(w, 2 * mu) + 40 * T + 2;

  const f = (eps) => {
    const dm = w - 2 * eps;
    if (Math.abs(dm) < 1e-7 * (w + 1)) return -Gp0 / (4 * w);
    return (grapheneG(eps, mu, T) - G0) / (dm * (w + 2 * eps));
  };

  let pts = [Math.max(0, mu - 10 * T), mu + 10 * T, 0.5 * w, w]
    .filter((x) => x > 0 && x < Lam);
  pts.push(0, Lam);
  pts.sort((a, b) => a - b);
  pts = pts.filter((x, i) => i === 0 || x - pts[i - 1] > 1e-12);

  let I = 0;
  for (let i = 0; i < pts.length - 1; i++) {
    const width = pts[i + 1] - pts[i];
    const n0 = Math.max(120, Math.ceil(40 * width / T));
    const nEven = n0 % 2 === 0 ? n0 : n0 + 1;
    I += simpson(f, pts[i], pts[i + 1], nEven);
  } 

  I += -((1 - G0) / (4 * w)) * Math.log((2 * Lam + w) / (2 * Lam - w));

  sRe += SIGMA0 * G0;
  sIm += SIGMA0 * (4 * w / Math.PI) * I;

  return { re: sRe, im: sIm };
}

function grapheneEffectiveIndex(sigma, wl_nm, d_nm) {
  const d = (d_nm || 0.34) * 1e-9;
  const omega = (2 * Math.PI * 2.99792458e8) / (wl_nm * 1e-9); // rad/s
  const s = 1 / (EPS0 * omega * d);
  // ε = 1 + iσs  →  {re: 1 − σ_im·s, im: σ_re·s}
  return cSqrt({ re: 1 - sigma.im * s, im: sigma.re * s });
}

// MARK: Non-Graphene TMM
// -----------------------------------------------------------------------------------------------

function cosThetaInLayer(nIncidentRe, theta0Rad, nLayer) {
  const beta = nIncidentRe * Math.sin(theta0Rad);
  const ratio = cDiv({ re: beta, im: 0 }, nLayer);
  const inside = cSub({ re: 1, im: 0 }, cMul(ratio, ratio));
  return cSqrt(inside);
}

function fresnelRT(ni, nj, cosI, cosJ, pol) {
  let r, t;
  if (pol === "s") {
    const num = cSub(cMul(ni, cosI), cMul(nj, cosJ));
    const den = cAdd(cMul(ni, cosI), cMul(nj, cosJ));
    r = cDiv(num, den);
    t = cDiv(cScaleR(cMul(ni, cosI), 2), den);
  } else {
    const num = cSub(cMul(nj, cosI), cMul(ni, cosJ));
    const den = cAdd(cMul(nj, cosI), cMul(ni, cosJ));
    r = cDiv(num, den);
    t = cDiv(cScaleR(cMul(ni, cosI), 2), den);
  }
  return [r, t];
}

function matMul2(A, B) {
  return [
    [cAdd(cMul(A[0][0], B[0][0]), cMul(A[0][1], B[1][0])), cAdd(cMul(A[0][0], B[0][1]), cMul(A[0][1], B[1][1]))],
    [cAdd(cMul(A[1][0], B[0][0]), cMul(A[1][1], B[1][0])), cAdd(cMul(A[1][0], B[0][1]), cMul(A[1][1], B[1][1]))],
  ];
}

// MARK: Graphene TMM
// -----------------------------------------------------------------------------------------------

function grapheneInterfaceMatrix(ni, nj, cosI, cosJ, pol, sigma) {
  const sEta = { re: sigma.re * ETA0, im: sigma.im * ETA0 };
  const one = { re: 1, im: 0 };
  if (pol === "s") {
    const q1 = cMul(ni, cosI);
    const eta = cDiv(cMul(nj, cosJ), q1);
    const xi = cDiv(sEta, q1);
    return [
      [cScaleR(cAdd(cAdd(one, eta), xi), 0.5), cScaleR(cAdd(cSub(one, eta), xi), 0.5)],
      [cScaleR(cSub(cSub(one, eta), xi), 0.5), cScaleR(cSub(cAdd(one, eta), xi), 0.5)],
    ];
  } else {
    const nu = cDiv(nj, ni);
    const gam = cDiv(cosJ, cosI);
    const zeta = cDiv(cMul(sEta, cosJ), ni);
    return [
      [cScaleR(cAdd(cAdd(nu, gam), zeta), 0.5), cScaleR(cSub(cSub(nu, gam), zeta), 0.5)],
      [cScaleR(cAdd(cSub(nu, gam), zeta), 0.5), cScaleR(cSub(cAdd(nu, gam), zeta), 0.5)],
    ];
  }
}

// MARK: Main Solver
// -----------------------------------------------------------------------------------------------

function solveTMM(wavelengths, nIncidentRe, nSubstrateRe, layers, theta0Deg, pol, grapheneParams) {
  const gp = Object.assign({}, GRAPHENE_DEFAULTS, grapheneParams || {});
  const T_eV = KB_EV * gp.T_K;

  const mediaLayers = [];
  const ifaceGraphene = [];
  const grapheneInfo = [];
  let pending = null;

  for (const L of layers) {
    if (L.type === "graphene") {
      const { n_m2, mu_eV } = grapheneGate(L.Vg || 0, gp);
      const { tau_s, Gamma_eV } = grapheneScattering(mu_eV, T_eV, gp);
      const desc = { mu_eV, Gamma_eV };
      desc.sigmaArr = wavelengths.map((wl) =>
        grapheneSigma(EV_NM / wl, mu_eV, T_eV, Gamma_eV)
      );
      grapheneInfo.push({
        Vg: L.Vg || 0,
        n_cm2: n_m2 * 1e-4,
        mu_eV,
        tau_fs: tau_s * 1e15,
        Gamma_meV: Gamma_eV * 1e3,
        lambdaPauli_nm: Math.abs(mu_eV) > 1e-4 ? EV_NM / (2 * Math.abs(mu_eV)) : Infinity,
      });
      pending = pending ? pending.concat([desc]) : [desc];
    } else {
      ifaceGraphene.push(pending);
      pending = null;
      mediaLayers.push(L);
    }
  }
  ifaceGraphene.push(pending);

  const theta0 = (theta0Deg * Math.PI) / 180;
  const nLambda = wavelengths.length;
  const R = new Array(nLambda);
  const T = new Array(nLambda);
  const A = new Array(nLambda);

  const nIncident = { re: nIncidentRe, im: 0 };
  const nSubstrate = { re: nSubstrateRe, im: 0 };

  for (let w = 0; w < nLambda; w++) {
    const wl = wavelengths[w];

    const nFull = [nIncident, ...mediaLayers.map((l) => l.n), nSubstrate];
    const cosList = nFull.map((n) => cosThetaInLayer(nIncidentRe, theta0, n));

    let M = [
      [{ re: 1, im: 0 }, { re: 0, im: 0 }],
      [{ re: 0, im: 0 }, { re: 1, im: 0 }],
    ];

    const nInterfaces = nFull.length - 1;
    for (let i = 0; i < nInterfaces; i++) {
      const ni = nFull[i];
      const nj = nFull[i + 1];
      const cosI = cosList[i];
      const cosJ = cosList[i + 1];

      let D;
      const sheets = ifaceGraphene[i];
      if (sheets && sheets.length) {
        let sig = { re: 0, im: 0 };
        for (const dsc of sheets) sig = cAdd(sig, dsc.sigmaArr[w]);
        D = grapheneInterfaceMatrix(ni, nj, cosI, cosJ, pol, sig);
      } else {
        const [r, t] = fresnelRT(ni, nj, cosI, cosJ, pol);
        D = [
          [cDiv({ re: 1, im: 0 }, t), cDiv(r, t)],
          [cDiv(r, t), cDiv({ re: 1, im: 0 }, t)],
        ];
      }
      M = matMul2(M, D);

      const isLastInterface = i === nInterfaces - 1;
      if (!isLastInterface) {
        const d = mediaLayers[i].d;
        const k0 = (2 * Math.PI) / wl;
        const phase = cScaleR(cMul(nj, cosJ), k0 * d);
        const iPhase = { re: -phase.im, im: phase.re };
        const negIPhase = { re: phase.im, im: -phase.re };
        const P = [
          [cExp(negIPhase), { re: 0, im: 0 }],
          [{ re: 0, im: 0 }, cExp(iPhase)],
        ];
        M = matMul2(M, P);
      }
    }

    const rTot = cDiv(M[1][0], M[0][0]);
    const tTot = cDiv({ re: 1, im: 0 }, M[0][0]);

    const Rw = cAbs2(rTot);
    const cos0 = cosList[0];
    const cosS = cosList[cosList.length - 1];
    const Tw = ((nSubstrateRe * cosS.re) / (nIncidentRe * cos0.re)) * cAbs2(tTot);

    R[w] = Rw;
    T[w] = Tw;
    A[w] = 1 - Rw - Tw;
  }

  return { R, T, A, graphene: grapheneInfo };
}

// MARK: Absorption Cross Checks
// -----------------------------------------------------------------------------------------------
function cInv2(M) {
  const det = cSub(cMul(M[0][0], M[1][1]), cMul(M[0][1], M[1][0]));
  return [
    [cDiv(M[1][1], det), cScaleR(cDiv(M[0][1], det), -1)],
    [cScaleR(cDiv(M[1][0], det), -1), cDiv(M[0][0], det)],
  ];
}

function computeFieldAbsorption(wavelengths, nIncidentRe, nSubstrateRe, layers, theta0Deg, pol, grapheneParams, sheetOccurrence) {
  if (pol !== 's') {
    throw new Error('computeFieldAbsorption: s-pol only (see p-pol basis note above).');
  }
  const occ = sheetOccurrence || 0;
  const gp = Object.assign({}, GRAPHENE_DEFAULTS, grapheneParams || {});
  const T_eV = KB_EV * gp.T_K;

  const mediaLayers = [];
  const ifaceGraphene = [];
  let pending = null;
  let sheetCount = -1;
  let targetInterfaceIdx = -1;
  let targetDescs = null;

  const flushInterface = () => {
    ifaceGraphene.push(pending);
    if (pending) {
      sheetCount++;
      if (sheetCount === occ) {
        targetInterfaceIdx = ifaceGraphene.length - 1;
        targetDescs = pending;
      }
    }
    pending = null;
  };

  for (const L of layers) {
    if (L.type === 'graphene') {
      const { mu_eV } = grapheneGate(L.Vg || 0, gp);
      const { Gamma_eV } = grapheneScattering(mu_eV, T_eV, gp);
      pending = pending ? pending.concat([{ mu_eV, Gamma_eV }]) : [{ mu_eV, Gamma_eV }];
    } else {
      flushInterface();
      mediaLayers.push(L);
    }
  }
  flushInterface();

  if (targetInterfaceIdx === -1) {
    throw new Error(`computeFieldAbsorption: no graphene sheet occurrence #${occ} found.`);
  }

  const theta0 = (theta0Deg * Math.PI) / 180;
  const nIncident = { re: nIncidentRe, im: 0 };
  const nSubstrate = { re: nSubstrateRe, im: 0 };
  const Afield = new Array(wavelengths.length);

  for (let w = 0; w < wavelengths.length; w++) {
    const wl = wavelengths[w];
    const nFull = [nIncident, ...mediaLayers.map(l => l.n), nSubstrate];
    const cosList = nFull.map(n => cosThetaInLayer(nIncidentRe, theta0, n));

    let M = [[{ re: 1, im: 0 }, { re: 0, im: 0 }], [{ re: 0, im: 0 }, { re: 1, im: 0 }]];
    let Mleft = null;
    const nInterfaces = nFull.length - 1;

    let sigTarget = { re: 0, im: 0 };
    for (const d of targetDescs) sigTarget = cAdd(sigTarget, grapheneSigma(EV_NM / wl, d.mu_eV, T_eV, d.Gamma_eV));

    for (let i = 0; i < nInterfaces; i++) {
      if (i === targetInterfaceIdx) Mleft = M; // snapshot BEFORE this interface's D

      const ni = nFull[i], nj = nFull[i + 1];
      const cosI = cosList[i], cosJ = cosList[i + 1];
      let D;
      const sheets = ifaceGraphene[i];
      if (sheets && sheets.length) {
        let sig = { re: 0, im: 0 };
        for (const dsc of sheets) sig = cAdd(sig, grapheneSigma(EV_NM / wl, dsc.mu_eV, T_eV, dsc.Gamma_eV));
        D = grapheneInterfaceMatrix(ni, nj, cosI, cosJ, pol, sig);
      } else {
        const [r, t] = fresnelRT(ni, nj, cosI, cosJ, pol);
        D = [[cDiv({ re: 1, im: 0 }, t), cDiv(r, t)], [cDiv(r, t), cDiv({ re: 1, im: 0 }, t)]];
      }
      M = matMul2(M, D);

      const isLast = i === nInterfaces - 1;
      if (!isLast) {
        const d = mediaLayers[i].d;
        const k0 = (2 * Math.PI) / wl;
        const phase = cScaleR(cMul(nj, cosJ), k0 * d);
        const iPhase = { re: -phase.im, im: phase.re };
        const negIPhase = { re: phase.im, im: -phase.re };
        const P = [[cExp(negIPhase), { re: 0, im: 0 }], [{ re: 0, im: 0 }, cExp(iPhase)]];
        M = matMul2(M, P);
      }
    }

    const rTot = cDiv(M[1][0], M[0][0]);
    const invMleft = cInv2(Mleft);
    const frontVec = [{ re: 1, im: 0 }, rTot];
    const Eplus = cAdd(cMul(invMleft[0][0], frontVec[0]), cMul(invMleft[0][1], frontVec[1]));
    const Eminus = cAdd(cMul(invMleft[1][0], frontVec[0]), cMul(invMleft[1][1], frontVec[1]));
    const Etan = cAdd(Eplus, Eminus);

    const cos0 = cosList[0];
    Afield[w] = (sigTarget.re * cAbs2(Etan) * ETA0) / (nIncidentRe * cos0.re);
  }

  return Afield;
}

function refinePeak(fn, guessWl, windowNm, tolNm) {
  let lo = guessWl - windowNm, hi = guessWl + windowNm;
  for (let iter = 0; iter < 40 && (hi - lo) > tolNm; iter++) {
    const m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
    if (fn(m1) < fn(m2)) lo = m1; else hi = m2;
  }
  return (lo + hi) / 2;
}

// MARK: console testing stuff
// -----------------------------------------------------------------------------------------------

/*
PRESETS[0].apply();
console.log(JSON.stringify(layerState, null, 1));
const solverLayers = layerState.map(l =>
  l.type === 'graphene' ? { type:'graphene', Vg: l.Vg||0 } : { n:{re:l.n, im:l.k||0}, d:l.d }
);
const critLayersVg0 = solverLayers;
const Aat = (wl) => solveTMM([wl], 1.0, 1.5, critLayersVg0, 0, 's').A[0];
const peakWl = refinePeak(Aat, 3000.2, 5, 1e-4);
console.log('refined peak wavelength =', peakWl, '  A =', Aat(peakWl));
const testWl = peakWl;
const AfieldAtPeak = computeFieldAbsorption([testWl], 1.0, 1.5, critLayersVg0, 0, 's');
const rtAtPeak = solveTMM([testWl], 1.0, 1.5, critLayersVg0, 0, 's').A[0];
console.log('A_field =', AfieldAtPeak[0], '   1-R-T =', rtAtPeak, '   diff =', Math.abs(AfieldAtPeak[0]-rtAtPeak));
*/

/*
function buildPhysicalCritLayers(lamC, nFront, nBack, Vg) {
  const dH = lamC/(4*N_HIGH), dL = lamC/(4*N_LOW);
  const mirror = (k) => Array.from({length:k}, (_,i) => i%2===0
    ? { n:{re:N_HIGH,im:0}, d:dH } : { n:{re:N_LOW,im:0}, d:dL });
  return [
    ...mirror(nFront),
    { n:{re:N_LOW,im:0}, d:dL/2 },   // physical: real material, half-width
    { type:'graphene', Vg },
    { n:{re:N_LOW,im:0}, d:dL/2 },
    ...mirror(nBack),
  ];
}
*/

/*
for (let nf = 6; nf <= 10; nf++) {
  for (let nb = 14; nb <= 18; nb++) {
    const layers = buildPhysicalCritLayers(3000, nf, nb, 0);
    const pk = refinePeak(wl => solveTMM([wl],1.0,1.5,layers,0,'s').A[0], 3000, 30, 1e-4);
    const A = solveTMM([pk],1.0,1.5,layers,0,'s').A[0];
    console.log(nf, nb, pk.toFixed(3), A.toFixed(5));
  }
}
*/