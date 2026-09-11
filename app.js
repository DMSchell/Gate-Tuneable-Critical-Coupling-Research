const CENTER_WL = 550.0;
const N_HIGH = 2.4, N_LOW = 1.45;
const D_HIGH = CENTER_WL / (4 * N_HIGH), D_LOW = CENTER_WL / (4 * N_LOW);

function defaultLayers() {
  const arr = [];
  for (let i = 0; i < 8; i++) {
    arr.push({ n: i % 2 === 0 ? N_HIGH : N_LOW, k: 0, d: i % 2 === 0 ? D_HIGH : D_LOW });
  }
  return arr;
}

let layerState = defaultLayers();
let result = { wavelengths: null, R: null, T: null, A: null, layers: [], graphene: [], warnings: [], gp: null, n_incident: 1.0, n_substrate: 1.5 };

// MARK: layer list UI
// ------------------------------------------------------------------------------------------

function renderLayers() {
  const container = document.getElementById('layers-container');
  container.innerHTML = '';

  layerState.forEach((layer, idx) => {
    const isG = layer.type === 'graphene';
    const row = document.createElement('div');
    row.className = 'layer-row';
    const typeSel = `
      <div><label>type</label>
        <select data-field="type" data-idx="${idx}">
          <option value="dielectric"${!isG ? ' selected' : ''}>diel.</option>
          <option value="graphene"${isG ? ' selected' : ''}>graph.</option>
        </select>
      </div>`;
    if (isG) {
      row.innerHTML = `
        <div class="layer-idx">${idx + 1}</div>
        ${typeSel}
        <div><label>V_g (V)</label><input type="number" step="1" value="${layer.Vg}" data-field="Vg" data-idx="${idx}"></div>
        <div class="graphene-note" style="grid-column: span 2;">monolayer, d = 0</div>
        <button class="remove-btn" data-idx="${idx}" title="Remove layer" aria-label="Remove layer ${idx + 1}">×</button>
      `;
    } else {
      row.innerHTML = `
        <div class="layer-idx">${idx + 1}</div>
        ${typeSel}
        <div><label>n</label><input type="number" step="0.01" value="${layer.n}" data-field="n" data-idx="${idx}"></div>
        <div><label>k (abs)</label><input type="number" step="0.001" value="${layer.k}" data-field="k" data-idx="${idx}"></div>
        <div><label>d (nm)</label><input type="number" step="1" value="${layer.d}" data-field="d" data-idx="${idx}"></div>
        <button class="remove-btn" data-idx="${idx}" title="Remove layer" aria-label="Remove layer ${idx + 1}">×</button>
      `;
    }
    container.appendChild(row);
  });

  container.querySelectorAll('input').forEach((inp) => {
    inp.addEventListener('change', (e) => {
      const idx = parseInt(e.target.dataset.idx);
      const field = e.target.dataset.field;
      layerState[idx][field] = parseFloat(e.target.value) || 0;
    });
  });

  container.querySelectorAll('select[data-field="type"]').forEach((sel) => {   // v2
    sel.addEventListener('change', (e) => {
      const idx = parseInt(e.target.dataset.idx);
      layerState[idx] = e.target.value === 'graphene'
        ? { type: 'graphene', Vg: 0 }
        : { n: 1.5, k: 0, d: 100 };
      renderLayers();
    });
  });

  container.querySelectorAll('.remove-btn').forEach((btn) => {
    btn.addEventListener('click', (e) => {
      const idx = parseInt(e.target.dataset.idx);
      layerState.splice(idx, 1);
      renderLayers();
    });
  });

  document.getElementById('layer-count').textContent = `(${layerState.length})`;
  populateVgLayerSelect();                                                     // v3
}

// MARK: Vg-sweep layer selector
// ------------------------------------------------------------------------------------------

function populateVgLayerSelect() {
  const sel = document.getElementById('vg-layer-select');
  if (!sel) return;
  const prevValue = sel.value;
  sel.innerHTML = '';
  layerState.forEach((l, idx) => {
    if (l.type === 'graphene') {
      const opt = document.createElement('option');
      opt.value = idx;
      opt.textContent = `Layer ${idx + 1} (currently V_g=${l.Vg} V)`;
      sel.appendChild(opt);
    }
  });
  if (sel.options.length === 0) {
    const opt = document.createElement('option');
    opt.value = '';
    opt.textContent = 'No graphene layer in the current stack';
    sel.appendChild(opt);
  } else if ([...sel.options].some((o) => o.value === prevValue)) {
    sel.value = prevValue; // keep selection stable across re-renders when possible
  }
}

function addLayer() {
  const last = layerState[layerState.length - 1];
  // v2: never auto-duplicate a graphene sheet (adjacent sheets are almost never intended)
  if (!last || last.type === 'graphene') layerState.push({ n: 1.5, k: 0, d: 100 });
  else layerState.push({ n: last.n, k: last.k, d: last.d });
  renderLayers();
  document.getElementById('layers-container').scrollTop = 1e6;
}

function resetToExample() {
  document.getElementById('n_incident').value = "1.0";
  document.getElementById('n_substrate').value = "1.5";
  document.getElementById('theta').value = "0";
  document.getElementById('pol').value = "s";
  document.getElementById('wl_min').value = "300";
  document.getElementById('wl_max').value = "900";
  document.getElementById('wl_pts').value = "400";
  layerState = defaultLayers();
  renderLayers();
  calculate();
}

// MARK: Presets
// ------------------------------------------------------------------------------------------

const PRESETS = [
  {
    id: 'main',
    name: 'Critical-coupling cavity (V_g = 32 V)',
    description: 'Asymmetric TiO₂/SiO₂ Bragg cavity (7 front / 21 back periods, λc = 3 µm) ' +
               'with a graphene monolayer at the antinode of a half-wave low-index defect. ' +
               'Gated to V_g = 32 V, at which δ = γ_e and simulated absorption reaches ' +
               'A ≈ 0.9987. See Table 1 and §4.1.',
    apply: () => {
      const lc = 3000;
      const dH = lc / (4 * N_HIGH), dL = lc / (4 * N_LOW);
      const mirror = (nPeriods) => {
        const m = [];
        for (let i = 0; i < nPeriods; i++) m.push(i % 2 === 0 ? { n: N_HIGH, k: 0, d: dH } : { n: N_LOW, k: 0, d: dL });
        return m;
      };
      layerState = [
        ...mirror(7),
        { n: N_LOW, k: 0, d: dL},
        { type: 'graphene', Vg: 32 },
        { n: N_LOW, k: 0, d: dL },
        ...mirror(21),
      ];
      document.getElementById('gr_TK').value = "300";
      document.getElementById('theta').value = "0";
      document.getElementById('wl_min').value = "2000";
      document.getElementById('wl_max').value = "4500";
      document.getElementById('wl_pts').value = "5000";
    },
  },
  {
    id: 'AtoT',
    name: 'Critical-coupling cavity (V_g = 32 V)',
    description: 'Asymmetric TiO₂/SiO₂ Bragg cavity (7 front / 21 back periods, λc = 3 µm) ' +
               'with a graphene monolayer at the antinode of a half-wave low-index defect. ' +
               'Gated to V_g = 32 V, at which δ = γ_e and simulated absorption reaches ' +
               'A ≈ 0.9987. See Table 1 and §4.1.',
    apply: () => {
      const lc = 3000;
      const dH = lc / (4 * N_HIGH), dL = lc / (4 * N_LOW);
      const mirror = (nPeriods) => {
        const m = [];
        for (let i = 0; i < nPeriods; i++) m.push(i % 2 === 0 ? { n: N_HIGH, k: 0, d: dH } : { n: N_LOW, k: 0, d: dL });
        return m;
      };
      layerState = [
        ...mirror(8),
        { n: N_LOW, k: 0, d: dL/2},
        { type: 'graphene', Vg: 32 },
        { n: N_LOW, k: 0, d: dL/2 },
        ...mirror(9),
      ];
      document.getElementById('gr_TK').value = "300";
      document.getElementById('theta').value = "0";
      document.getElementById('wl_min').value = "2000";
      document.getElementById('wl_max').value = "4500";
      document.getElementById('wl_pts').value = "5000";
    },
  },
  {
    id: 'freestanding',
    name: 'Freestanding Graphene',
    description: 'Freestanding graphene sheet (n = 1.0, d = 0) in air. See UNKNOWN.',
    apply: () => {
      layerState = [
        { n: 1.0, k: 0, d: 30},
        { type: 'graphene', Vg: 40 },
        { n: 1.0, k: 0, d: 30 },
      ];
      document.getElementById('n_incident').value = "1.0";
      document.getElementById('n_substrate').value = "1.0";
      document.getElementById('gr_TK').value = "300";
      document.getElementById('theta').value = "0";
      document.getElementById('wl_min').value = "2000";
      document.getElementById('wl_max').value = "4500";
      document.getElementById('wl_pts').value = "5000";
    },
  },
];

function applyPreset(preset) {
  preset.apply();
  renderLayers();
  calculate();
  closePresetMenu();
}

function renderPresetMenu() {
  const menu = document.getElementById('preset-menu');
  menu.innerHTML = '';
  PRESETS.forEach((p) => {
    const btn = document.createElement('button');
    btn.className = 'preset-item';
    btn.type = 'button';
    btn.innerHTML = `<span class="preset-name">${p.name}</span><span class="preset-desc">${p.description}</span>`;
    btn.addEventListener('click', () => applyPreset(p));
    menu.appendChild(btn);
  });
}

function openPresetMenu() {
  const menu = document.getElementById('preset-menu');
  const btn = document.getElementById('preset-btn');
  menu.hidden = false;
  btn.setAttribute('aria-expanded', 'true');
  document.addEventListener('click', onDocClickForPresetMenu, true);
  document.addEventListener('keydown', onKeydownForPresetMenu, true);
}

function closePresetMenu() {
  const menu = document.getElementById('preset-menu');
  const btn = document.getElementById('preset-btn');
  menu.hidden = true;
  btn.setAttribute('aria-expanded', 'false');
  document.removeEventListener('click', onDocClickForPresetMenu, true);
  document.removeEventListener('keydown', onKeydownForPresetMenu, true);
}

function onDocClickForPresetMenu(e) {
  const dropdown = document.getElementById('preset-dropdown');
  if (!dropdown.contains(e.target)) closePresetMenu();
}

function onKeydownForPresetMenu(e) {
  if (e.key === 'Escape') closePresetMenu();
}

document.getElementById('preset-btn').addEventListener('click', (e) => {
  e.stopPropagation();
  const menu = document.getElementById('preset-menu');
  if (menu.hidden) openPresetMenu(); else closePresetMenu();
});

renderPresetMenu();

// MARK: Small Helpers
// ------------------------------------------------------------------------------------------

function linspace(a, b, n) {
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + (b - a) * i / (n - 1);
  return arr;
}

function nearestIndex(arr, v) {
  let lo = 0, hi = arr.length - 1;
  if (v <= arr[0]) return 0;
  if (v >= arr[hi]) return hi;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] < v) lo = mid; else hi = mid;
  }
  return (v - arr[lo] < arr[hi] - v) ? lo : hi;
}

// v2: read graphene/gate inputs (isFinite, not ||, so 0 is a legal value)
function collectGrapheneParams() {
  const v = (id, dflt) => {
    const x = parseFloat(document.getElementById(id).value);
    return isFinite(x) ? x : dflt;
  };
  return {
    epsOx: v('gr_epsOx', 3.9), tox_nm: v('gr_tox', 285), VDirac: v('gr_vdirac', 0),
    T_K: v('gr_TK', 300), nimp_cm2: v('gr_nimp', 1e11), epsG: v('gr_epsG', 2.45),
    vF: v('gr_vF', 1.0e6),
  };
}

function graphenePlacementWarnings() {
  const w = [];
  layerState.forEach((l, i) => {
    if (l.type !== 'graphene') return;
    if (layerState[i + 1] && layerState[i + 1].type === 'graphene')
      w.push(`Layers ${i + 1}–${i + 2}: adjacent graphene sheets are coincident (their conductivities add).`);
  });
  return w;
}

// MARK: Spectrum Chart
// ------------------------------------------------------------------------------------------

function drawChart() {
  const svg = document.getElementById('chart');
  svg.innerHTML = '';
  const W = 720, H = 340, padL = 46, padR = 14, padT = 14, padB = 34;
  const plotW = W - padL - padR, plotH = H - padT - padB;
  const ns = 'http://www.w3.org/2000/svg';

  if (!result.wavelengths) {
    const t = document.createElementNS(ns, 'text');
    t.setAttribute('x', W / 2); t.setAttribute('y', H / 2);
    t.setAttribute('text-anchor', 'middle'); t.setAttribute('fill', '#6b7280');
    t.textContent = 'Press Calculate to run the simulation';
    svg.appendChild(t);
    return;
  }

  const wl = result.wavelengths;
  const wlMin = wl[0], wlMax = wl[wl.length - 1];
  const xScale = (v) => padL + (v - wlMin) / (wlMax - wlMin) * plotW;
  const yScale = (v) => padT + (1 - v) * plotH;

  for (let f = 0; f <= 1.0001; f += 0.2) {
    const y = yScale(f);
    const line = document.createElementNS(ns, 'line');
    line.setAttribute('x1', padL); line.setAttribute('x2', W - padR);
    line.setAttribute('y1', y); line.setAttribute('y2', y);
    line.setAttribute('stroke', '#e5e8ec'); line.setAttribute('stroke-width', '1');
    svg.appendChild(line);

    const lbl = document.createElementNS(ns, 'text');
    lbl.setAttribute('x', padL - 8); lbl.setAttribute('y', y + 4);
    lbl.setAttribute('text-anchor', 'end'); lbl.setAttribute('font-size', '11'); lbl.setAttribute('fill', '#6b7280');
    lbl.textContent = f.toFixed(1);
    svg.appendChild(lbl);
  }

  for (let i = 0; i <= 4; i++) {
    const v = wlMin + (wlMax - wlMin) * i / 4;
    const x = xScale(v);
    const lbl = document.createElementNS(ns, 'text');
    lbl.setAttribute('x', x); lbl.setAttribute('y', H - padB + 18);
    lbl.setAttribute('text-anchor', 'middle'); lbl.setAttribute('font-size', '11'); lbl.setAttribute('fill', '#6b7280');
    lbl.textContent = Math.round(v);
    svg.appendChild(lbl);
  }

  const xlabel = document.createElementNS(ns, 'text');
  xlabel.setAttribute('x', W / 2); xlabel.setAttribute('y', H - 4);
  xlabel.setAttribute('text-anchor', 'middle'); xlabel.setAttribute('font-size', '12'); xlabel.setAttribute('fill', '#6b7280');
  xlabel.textContent = 'Wavelength (nm)';
  svg.appendChild(xlabel);

  function pathFor(arr) {
    let d = '';
    for (let i = 0; i < wl.length; i++) {
      const x = xScale(wl[i]), y = yScale(Math.max(0, Math.min(1, arr[i])));
      d += (i === 0 ? 'M' : 'L') + x.toFixed(2) + ',' + y.toFixed(2) + ' ';
    }
    return d;
  }

  const css = getComputedStyle(document.documentElement);

  const pathA = document.createElementNS(ns, 'path');
  pathA.setAttribute('d', pathFor(result.A));
  pathA.setAttribute('fill', 'none'); pathA.setAttribute('stroke', css.getPropertyValue('--accent-a'));
  pathA.setAttribute('stroke-width', '1.5'); pathA.setAttribute('stroke-dasharray', '4 3');
  svg.appendChild(pathA);

  const pathT = document.createElementNS(ns, 'path');
  pathT.setAttribute('d', pathFor(result.T));
  pathT.setAttribute('fill', 'none'); pathT.setAttribute('stroke', css.getPropertyValue('--accent-t'));
  pathT.setAttribute('stroke-width', '2');
  svg.appendChild(pathT);

  const pathR = document.createElementNS(ns, 'path');
  pathR.setAttribute('d', pathFor(result.R));
  pathR.setAttribute('fill', 'none'); pathR.setAttribute('stroke', css.getPropertyValue('--accent-r'));
  pathR.setAttribute('stroke-width', '2');
  svg.appendChild(pathR);

  const probeWl = parseFloat(document.getElementById('probe-slider').value);
  const px = xScale(probeWl);
  const probeLine = document.createElementNS(ns, 'line');
  probeLine.setAttribute('x1', px); probeLine.setAttribute('x2', px);
  probeLine.setAttribute('y1', padT); probeLine.setAttribute('y2', H - padB);
  probeLine.setAttribute('stroke', '#1c2430'); probeLine.setAttribute('stroke-width', '1');
  probeLine.setAttribute('stroke-dasharray', '2 2');
  svg.appendChild(probeLine);
}

drawChart();

// MARK: Vg sweep
// ------------------------------------------------------------------------------------------
// ---- v3: Vg sweep, directly on the live layer stack --------------------
// Sweeps V_g on ONE selected graphene layer, holding every other layer
// (mirrors, other graphene sheets, global n_incident/n_substrate/theta/pol)
// exactly as currently configured in the UI. Reuses the SAME conversion
// calculate() applies to layerState — it's your actual configured stack,
// just looped over Vg. (Replaces the earlier buildCavity()-based sweep,
// which rebuilt its own separate geometry from N_front/N_back rather than
// using the live layer list.)

let vgSweepData = [];

function runLayerVgSweep() {
  const sel = document.getElementById('vg-layer-select');
  const targetIdx = parseInt(sel.value);
  const statusEl = document.getElementById('vg-status');
  const runBtn = document.getElementById('vg-run-btn');

  if (!Number.isFinite(targetIdx)) {
    statusEl.textContent = 'No graphene layer in the current stack — add one first.';
    return;
  }

  const n_incident = parseFloat(document.getElementById('n_incident').value) || 1.0;
  const n_substrate = parseFloat(document.getElementById('n_substrate').value) || 1.5;
  const theta = parseFloat(document.getElementById('theta').value) || 0;
  const pol = document.getElementById('pol').value;
  const gp = collectGrapheneParams();

  const wlProbe = parseFloat(document.getElementById('vg_wl_probe').value) || 550;
  const vgMin = parseFloat(document.getElementById('vg_min').value) || 0;
  const vgMax = parseFloat(document.getElementById('vg_max').value) || 100;
  const nPtsSweep = Math.max(parseInt(document.getElementById('vg_pts').value) || 50, 2);
  const relocate = document.getElementById('vg_relocate').value === 'yes';
  const windowNm = parseFloat(document.getElementById('vg_window').value) || 30;

  runBtn.disabled = true;
  statusEl.textContent = 'Running…';
  vgSweepData = [];
  let lastGuess = wlProbe;

  setTimeout(() => {
    for (let i = 0; i < nPtsSweep; i++) {
      const Vg = vgMin + (vgMax - vgMin) * i / (nPtsSweep - 1);

      const solverLayers = layerState.map((l, idx) => {
        if (l.type === 'graphene') {
          return { type: 'graphene', Vg: idx === targetIdx ? Vg : (l.Vg || 0) };
        }
        return { n: { re: l.n, im: l.k || 0 }, d: l.d };
      });

      let point;
      if (relocate) {
        const wl = linspace(lastGuess - windowNm, lastGuess + windowNm, 121);
        const r = solveTMM(wl, n_incident, n_substrate, solverLayers, theta, pol, gp);
        let best = 0;
        for (let k = 1; k < r.A.length; k++) if (r.A[k] > r.A[best]) best = k;
        point = { wl: wl[best], R: r.R[best], T: r.T[best], A: r.A[best] };
        lastGuess = point.wl;
      } else {
        const r = solveTMM([wlProbe], n_incident, n_substrate, solverLayers, theta, pol, gp);
        point = { wl: wlProbe, R: r.R[0], T: r.T[0], A: r.A[0] };
      }

      const { mu_eV } = grapheneGate(Vg, gp);
      vgSweepData.push({ Vg, wl: point.wl, R: point.R, T: point.T, A: point.A, mu_eV });
    }

    statusEl.textContent = `Done — ${nPtsSweep} points, V_g ${vgMin}–${vgMax} V, sweeping layer ${targetIdx + 1}.`;
    runBtn.disabled = false;

    const slider = document.getElementById('vg-probe-slider');
    slider.min = 0; slider.max = vgSweepData.length - 1; slider.value = 0;
    slider.disabled = false;

    drawVgChart();
    updateVgReadout(0);
  }, 20);
}

function drawVgChart() {
  const svg = document.getElementById('vg-chart');
  svg.innerHTML = '';
  const ns = 'http://www.w3.org/2000/svg';
  const W = 720, H = 300, padL = 46, padR = 14, padT = 14, padB = 34;
  const plotW = W - padL - padR, plotH = H - padT - padB;

  if (!vgSweepData.length) {
    const t = document.createElementNS(ns, 'text');
    t.setAttribute('x', W / 2); t.setAttribute('y', H / 2);
    t.setAttribute('text-anchor', 'middle'); t.setAttribute('fill', '#6b7280');
    t.textContent = 'Run a V_g sweep to see R/T/A';
    svg.appendChild(t);
    return;
  }

  const vgs = vgSweepData.map((p) => p.Vg);
  const vgMin = vgs[0], vgMax = vgs[vgs.length - 1];
  const xScale = (v) => padL + (v - vgMin) / ((vgMax - vgMin) || 1) * plotW;
  const yScale = (v) => padT + (1 - v) * plotH;

  for (let f = 0; f <= 1.0001; f += 0.2) {
    const y = yScale(f);
    const line = document.createElementNS(ns, 'line');
    line.setAttribute('x1', padL); line.setAttribute('x2', W - padR);
    line.setAttribute('y1', y); line.setAttribute('y2', y);
    line.setAttribute('stroke', '#e5e8ec'); line.setAttribute('stroke-width', '1');
    svg.appendChild(line);
    const lbl = document.createElementNS(ns, 'text');
    lbl.setAttribute('x', padL - 8); lbl.setAttribute('y', y + 4);
    lbl.setAttribute('text-anchor', 'end'); lbl.setAttribute('font-size', '11'); lbl.setAttribute('fill', '#6b7280');
    lbl.textContent = f.toFixed(1);
    svg.appendChild(lbl);
  }
  for (let i = 0; i <= 5; i++) {
    const v = vgMin + (vgMax - vgMin) * i / 5;
    const x = xScale(v);
    const lbl = document.createElementNS(ns, 'text');
    lbl.setAttribute('x', x); lbl.setAttribute('y', H - padB + 18);
    lbl.setAttribute('text-anchor', 'middle'); lbl.setAttribute('font-size', '11'); lbl.setAttribute('fill', '#6b7280');
    lbl.textContent = v.toFixed(1);
    svg.appendChild(lbl);
  }
  const xlabel = document.createElementNS(ns, 'text');
  xlabel.setAttribute('x', W / 2); xlabel.setAttribute('y', H - 4);
  xlabel.setAttribute('text-anchor', 'middle'); xlabel.setAttribute('font-size', '12'); xlabel.setAttribute('fill', '#6b7280');
  xlabel.textContent = 'V_g (V)';
  svg.appendChild(xlabel);

  function pathFor(key) {
    let d = '';
    vgSweepData.forEach((p, i) => {
      const x = xScale(p.Vg), y = yScale(Math.max(0, Math.min(1, p[key])));
      d += (i === 0 ? 'M' : 'L') + x.toFixed(2) + ',' + y.toFixed(2) + ' ';
    });
    return d;
  }

  const css = getComputedStyle(document.documentElement);

  const pathA = document.createElementNS(ns, 'path');
  pathA.setAttribute('d', pathFor('A'));
  pathA.setAttribute('fill', 'none'); pathA.setAttribute('stroke', css.getPropertyValue('--accent-a'));
  pathA.setAttribute('stroke-width', '2'); pathA.setAttribute('stroke-dasharray', '4 3');
  svg.appendChild(pathA);

  const pathT = document.createElementNS(ns, 'path');
  pathT.setAttribute('d', pathFor('T'));
  pathT.setAttribute('fill', 'none'); pathT.setAttribute('stroke', css.getPropertyValue('--accent-t'));
  pathT.setAttribute('stroke-width', '2');
  svg.appendChild(pathT);

  const pathR = document.createElementNS(ns, 'path');
  pathR.setAttribute('d', pathFor('R'));
  pathR.setAttribute('fill', 'none'); pathR.setAttribute('stroke', css.getPropertyValue('--accent-r'));
  pathR.setAttribute('stroke-width', '2');
  svg.appendChild(pathR);

  const idx = parseInt(document.getElementById('vg-probe-slider').value) || 0;
  const px = xScale(vgSweepData[idx].Vg);
  const probeLine = document.createElementNS(ns, 'line');
  probeLine.setAttribute('x1', px); probeLine.setAttribute('x2', px);
  probeLine.setAttribute('y1', padT); probeLine.setAttribute('y2', H - padB);
  probeLine.setAttribute('stroke', '#1c2430'); probeLine.setAttribute('stroke-width', '1');
  probeLine.setAttribute('stroke-dasharray', '2 2');
  svg.appendChild(probeLine);
}

function updateVgReadout(idx) {
  if (!vgSweepData.length) return;
  const p = vgSweepData[idx];
  document.getElementById('vg-probe-value').textContent =
    `V_g=${p.Vg.toFixed(2)} V, λ=${p.wl.toFixed(2)} nm, R=${p.R.toFixed(4)} T=${p.T.toFixed(4)} A=${p.A.toFixed(4)}`;

  let maxIdx = 0, minIdx = 0;
  for (let i = 1; i < vgSweepData.length; i++) {
    if (vgSweepData[i].A > vgSweepData[maxIdx].A) maxIdx = i;
    if (vgSweepData[i].A < vgSweepData[minIdx].A) minIdx = i;
  }
  const on = vgSweepData[maxIdx], off = vgSweepData[minIdx];
  const contrast = on.A / Math.max(off.A, 1e-9);

  const html =
    `Peak A = ${on.A.toFixed(4)} at V_g=${on.Vg.toFixed(2)} V (μ_c=${on.mu_eV.toFixed(4)} eV, λ=${on.wl.toFixed(1)} nm)<br>` +
    `Min A &nbsp;= ${off.A.toFixed(4)} at V_g=${off.Vg.toFixed(2)} V (μ_c=${off.mu_eV.toFixed(4)} eV, λ=${off.wl.toFixed(1)} nm)<br>` +
    `On/off contrast (max A / min A): ${contrast.toFixed(2)}×`;
  document.getElementById('vg-readout').innerHTML = html;
}

// MARK: beam diagram
// ------------------------------------------------------------------------------------------

function drawDiagram() {
  const svg = document.getElementById('diagram');
  svg.innerHTML = '';
  const ns = 'http://www.w3.org/2000/svg';
  const W = 720, H = 160;

  if (!result.wavelengths || result.layers.length === 0) {
    const t = document.createElementNS(ns, 'text');
    t.setAttribute('x', W / 2); t.setAttribute('y', H / 2);
    t.setAttribute('text-anchor', 'middle'); t.setAttribute('fill', '#6b7280');
    t.textContent = 'Press Calculate to build the diagram';
    svg.appendChild(t);
    return;
  }

  const probeWl = parseFloat(document.getElementById('probe-slider').value);
  const idx = nearestIndex(result.wavelengths, probeWl);
  const Rh = result.R[idx], Th = result.T[idx], Ah = result.A[idx];

  const xStart = 140, xEnd = 580, yTop = 30, yH = 70;
  const layers = result.layers;

  // v2: widths computed over dielectric layers only; graphene has zero width
  const dielectrics = layers.filter((l) => l.type !== 'graphene');
  const thicknesses = dielectrics.map((l) => l.d);
  let totalD = thicknesses.reduce((a, b) => a + b, 0) || 1;
  const minFrac = 0.03;
  let fracs = thicknesses.map((d) => Math.max(d / totalD, minFrac));
  const fracSum = fracs.reduce((a, b) => a + b, 0) || 1;
  fracs = fracs.map((f) => f / fracSum);

  const nReals = [result.n_incident, ...dielectrics.map((l) => l.n), result.n_substrate];
  const nMin = Math.min(...nReals), nMax = Math.max(...nReals);
  const span = Math.max(nMax - nMin, 1e-6);
  function colorFor(n) {
    const t = (n - nMin) / span;
    const lo = [253, 219, 199], hi = [140, 60, 10];
    const c = lo.map((v, i) => Math.round(v + (hi[i] - v) * t));
    return `rgb(${c[0]},${c[1]},${c[2]})`;
  }

  const incidentRect = document.createElementNS(ns, 'rect');
  incidentRect.setAttribute('x', 0); incidentRect.setAttribute('y', yTop);
  incidentRect.setAttribute('width', xStart); incidentRect.setAttribute('height', yH);
  incidentRect.setAttribute('fill', '#eef0f3'); incidentRect.setAttribute('stroke', '#c7ccd3');
  svg.appendChild(incidentRect);

  const incidentLbl = document.createElementNS(ns, 'text');
  incidentLbl.setAttribute('x', xStart / 2); incidentLbl.setAttribute('y', yTop + yH / 2 + 4);
  incidentLbl.setAttribute('text-anchor', 'middle'); incidentLbl.setAttribute('font-size', '11'); incidentLbl.setAttribute('fill', '#3a3f47');
  incidentLbl.textContent = `n0=${result.n_incident.toFixed(2)}`;
  svg.appendChild(incidentLbl);

  const subRect = document.createElementNS(ns, 'rect');
  subRect.setAttribute('x', xEnd); subRect.setAttribute('y', yTop);
  subRect.setAttribute('width', W - xEnd); subRect.setAttribute('height', yH);
  subRect.setAttribute('fill', '#eef0f3'); subRect.setAttribute('stroke', '#c7ccd3');
  svg.appendChild(subRect);

  const subLbl = document.createElementNS(ns, 'text');
  subLbl.setAttribute('x', xEnd + (W - xEnd) / 2); subLbl.setAttribute('y', yTop + yH / 2 + 4);
  subLbl.setAttribute('text-anchor', 'middle'); subLbl.setAttribute('font-size', '11'); subLbl.setAttribute('fill', '#3a3f47');
  subLbl.textContent = `ns=${result.n_substrate.toFixed(2)}`;
  svg.appendChild(subLbl);

  let x = xStart;
  const stackWidth = xEnd - xStart;
  let di = 0;
  layers.forEach((layer) => {
    if (layer.type === 'graphene') {
      const gl = document.createElementNS(ns, 'line');
      gl.setAttribute('x1', x); gl.setAttribute('x2', x);
      gl.setAttribute('y1', yTop - 6); gl.setAttribute('y2', yTop + yH + 6);
      gl.setAttribute('stroke', '#111'); gl.setAttribute('stroke-width', '2.5');
      svg.appendChild(gl);
      const gt = document.createElementNS(ns, 'text');
      gt.setAttribute('x', x); gt.setAttribute('y', yTop - 10);
      gt.setAttribute('text-anchor', 'middle'); gt.setAttribute('font-size', '10'); gt.setAttribute('fill', '#111');
      gt.textContent = `G (${layer.Vg} V)`;
      svg.appendChild(gt);
      return;
    }
    const w = fracs[di++] * stackWidth;
    const rect = document.createElementNS(ns, 'rect');
    rect.setAttribute('x', x); rect.setAttribute('y', yTop);
    rect.setAttribute('width', w); rect.setAttribute('height', yH);
    rect.setAttribute('fill', colorFor(layer.n)); rect.setAttribute('stroke', '#8a5a2a'); rect.setAttribute('stroke-width', '0.5');
    svg.appendChild(rect);

    if (w > 18) {
      const lbl = document.createElementNS(ns, 'text');
      lbl.setAttribute('x', x + w / 2); lbl.setAttribute('y', yTop + yH / 2 + 4);
      lbl.setAttribute('text-anchor', 'middle'); lbl.setAttribute('font-size', '10'); lbl.setAttribute('fill', '#2c1c08');
      lbl.textContent = layer.n.toFixed(2);
      svg.appendChild(lbl);
    }
    x += w;
  });

  function arrow(x1, y1, x2, y2, color, opacity) {
    const g = document.createElementNS(ns, 'g');
    g.setAttribute('opacity', opacity);

    const line = document.createElementNS(ns, 'line');
    line.setAttribute('x1', x1); line.setAttribute('y1', y1);
    line.setAttribute('x2', x2); line.setAttribute('y2', y2);
    line.setAttribute('stroke', color); line.setAttribute('stroke-width', '3');
    g.appendChild(line);

    const angle = Math.atan2(y2 - y1, x2 - x1);
    const headLen = 9;
    const p1x = x2 - headLen * Math.cos(angle - 0.4), p1y = y2 - headLen * Math.sin(angle - 0.4);
    const p2x = x2 - headLen * Math.cos(angle + 0.4), p2y = y2 - headLen * Math.sin(angle + 0.4);
    const head = document.createElementNS(ns, 'polygon');
    head.setAttribute('points', `${x2},${y2} ${p1x},${p1y} ${p2x},${p2y}`);
    head.setAttribute('fill', color);
    g.appendChild(head);

    svg.appendChild(g);
  }

  arrow(10, yTop + yH / 1.5, xStart - 8, yTop + yH / 1.5, '#1c2430', 1.0);
  arrow(xStart - 8, yTop + yH / 3, 10, yTop + yH / 3, '#c0392b', 0.15 + 0.85 * Rh);
  arrow(xEnd + 8, yTop + yH / 1.5, W - 10, yTop + yH / 1.5, '#2166ac', 0.15 + 0.85 * Th);

  const info = document.createElementNS(ns, 'text');
  info.setAttribute('x', W / 2); info.setAttribute('y', 10);
  info.setAttribute('text-anchor', 'middle'); info.setAttribute('font-size', '12'); info.setAttribute('font-weight', '600'); info.setAttribute('fill', '#1c2430');
  info.textContent = `λ = ${result.wavelengths[idx].toFixed(1)} nm    R = ${Rh.toFixed(3)}    T = ${Th.toFixed(3)}    A = ${Ah.toFixed(3)}`;
  svg.appendChild(info);
}

// MARK: text readout
// ------------------------------------------------------------------------------------------

function updateReadout() {
  if (!result.wavelengths) return;
  const probeWl = parseFloat(document.getElementById('probe-slider').value);
  const idx = nearestIndex(result.wavelengths, probeWl);

  let RPeakIdx = 0, TPeakIdx = 0, APeakIdx = 0;
  for (let i = 1; i < result.R.length; i++) if (result.R[i] > result.R[RPeakIdx]) RPeakIdx = i;
  for (let i = 1; i < result.T.length; i++) if (result.T[i] > result.T[TPeakIdx]) TPeakIdx = i;
  for (let i = 1; i < result.A.length; i++) if (result.A[i] > result.A[APeakIdx]) APeakIdx = i;

  let html =
    `Probe: λ = ${result.wavelengths[idx].toFixed(1)} nm &nbsp; R = ${result.R[idx].toFixed(4)} &nbsp; T = ${result.T[idx].toFixed(4)} &nbsp; A = ${result.A[idx].toFixed(4)}<br>` +
    `Peak reflectance: R_max = ${result.R[RPeakIdx].toFixed(4)} at λ = ${result.wavelengths[RPeakIdx].toFixed(1)} nm <br>` + 
    `Peak transmission: T_max = ${result.T[TPeakIdx].toFixed(4)} at λ = ${result.wavelengths[TPeakIdx].toFixed(1)} nm <br>` +
    `Peak absorption: A_max = ${result.A[APeakIdx].toFixed(4)} at λ = ${result.wavelengths[APeakIdx].toFixed(1)} nm <br>`;

  if (result.graphene && result.graphene.length) {                             // v2
    const T_eV = KB_EV * ((result.gp && result.gp.T_K) || 300);
    const probeE = EV_NM / result.wavelengths[idx];
    result.graphene.forEach((g, k) => {
      const s = grapheneSigma(probeE, g.mu_eV, T_eV, g.Gamma_meV * 1e-3);
      const sr = (s.re / SIGMA0), si = (s.im / SIGMA0);
      const lam = isFinite(g.lambdaPauli_nm) ? `${(g.lambdaPauli_nm / 1000).toFixed(2)} µm` : '—';
      html += `<br>G${k + 1}: V_g=${g.Vg} V | n=${g.n_cm2.toExponential(2)} cm⁻² | E_F=${g.mu_eV.toFixed(3)} eV | ` +
              `2E_F edge (P-B point): ${lam} | τ=${g.tau_fs.toFixed(0)} fs | σ/σ₀ @probe = ${sr.toFixed(2)} ${si < 0 ? '−' : '+'} ${Math.abs(si).toFixed(2)}i | ` +
              `Cg = ${(((8.8541878128e-12 * result.gp.epsOx) / (result.gp.tox_nm * 1e-9)) * 1e6).toFixed(1)} aF/µm² | µ = ${g.mu_eV.toFixed(3)} eV` ;
    });
  }
  if (result.warnings && result.warnings.length) {                             // v2
    html += `<br><span class="warn">⚠ ${result.warnings.join('<br>⚠ ')}</span>`;
  }

  document.getElementById('readout').innerHTML = html;
}

// ---- Calculate ----

function calculate() {
  console.log("Calculating...");
  console.time("Calculation-timer");
  const n_incident = parseFloat(document.getElementById('n_incident').value) || 1.0;
  const n_substrate = parseFloat(document.getElementById('n_substrate').value) || 1.5;
  const theta = parseFloat(document.getElementById('theta').value) || 0;
  const pol = document.getElementById('pol').value;
  const wlMin = parseFloat(document.getElementById('wl_min').value) || 300;
  const wlMax = parseFloat(document.getElementById('wl_max').value) || 900;
  const wlPts = Math.max(parseInt(document.getElementById('wl_pts').value) || 400, 2);

  // v2: two layer kinds go to the solver
  const solverLayers = layerState.map((l) =>
    l.type === 'graphene'
      ? { type: 'graphene', Vg: l.Vg || 0 }
      : { n: { re: l.n, im: l.k || 0 }, d: l.d }
  );
  const gp = collectGrapheneParams();
  const wavelengths = linspace(wlMin, wlMax, wlPts);

  const out = solveTMM(wavelengths, n_incident, n_substrate, solverLayers, theta, pol, gp);

  result = {
    wavelengths, R: out.R, T: out.T, A: out.A,
    graphene: out.graphene || [],
    warnings: graphenePlacementWarnings(),
    gp,
    layers: layerState.map((l) => (l.type === 'graphene' ? { type: 'graphene', Vg: l.Vg || 0 } : { n: l.n, k: l.k, d: l.d })),
    n_incident, n_substrate,
  };

  const slider = document.getElementById('probe-slider');
  slider.min = wlMin; slider.max = wlMax;
  slider.value = (wlMin + wlMax) / 2;
  document.getElementById('probe-value').textContent = `${Math.round(slider.value)} nm`;

  drawChart();
  drawDiagram();
  updateReadout();
  console.timeEnd("Calculation-timer");
}


// MARK: validation
// ------------------------------------------------------------------------------------------
function runGrapheneValidation() {
  const out = {};
  const air = { re: 1, im: 0 };

  { // 1. regression: no graphene → legacy path must match v1
    const wl = linspace(300, 900, 400);
    const r = solveTMM(wl, 1.0, 1.5, defaultLayers().map((l) => ({ n: { re: l.n, im: l.k }, d: l.d })), 0, 's');
    out['1. Bragg R(550) [≈ v1 value]'] = r.R[nearestIndex(wl, 550)].toFixed(4);
  }
  { // 2. freestanding sheet, Vg=0, θ=0
    const L = [{ n: air, d: 100 }, { type: 'graphene', Vg: 0 }, { n: air, d: 100 }];
    const idealParams = { T_K: 1e-3, nimp_cm2: 1e-30 };
    const s = solveTMM([550], 1, 1, L, 0, 's', idealParams);
    const p = solveTMM([550], 1, 1, L, 0, 'p', idealParams);
    out['2a. freestanding A [0.02241]'] = s.A[0].toFixed(5);
    out['2b. |A_s − A_p| [≈0]'] = Math.abs(s.A[0] - p.A[0]).toExponential(2);
  }
  { // 3. Kubo vs T→0 closed form (μ=0.3 eV, ħω=1 eV) — catches the factor-of-4
    const sg = grapheneSigma(1.0, 0.3, 1e-5, 1e-9);
    const imExpect = 4 * 0.3 / Math.PI - Math.log(1.6 / 0.4) / Math.PI;
    out['3a. Re σ/σ0 [1.000]'] = (sg.re / SIGMA0).toFixed(4);
    out[`3b. Im σ/σ0 [${imExpect.toFixed(4)}]`] = (sg.im / SIGMA0).toFixed(4);
  }
  { // 4. boundary-condition vs 0.34 nm effective-medium film, inside the cavity
    const gp = GRAPHENE_DEFAULTS, T_eV = KB_EV * gp.T_K, Vg = 60;
    const { mu_eV } = grapheneGate(Vg, gp);
    const { Gamma_eV } = grapheneScattering(mu_eV, T_eV, gp);
    const H = { re: N_HIGH, im: 0 }, Lo = { re: N_LOW, im: 0 };
    const dH = 3000 / (4 * N_HIGH), dL = 3000 / (4 * N_LOW);
    const mirror = [];
    for (let i = 0; i < 7; i++) mirror.push({ n: i % 2 ? Lo : H, d: i % 2 ? dL : dH });
    let maxDiff = 0;
    for (let k = 0; k <= 160; k++) {
      const wl = 2600 + k * 5;
      const bc = solveTMM([wl], 1, 1.5,
        [...mirror, { n: Lo, d: dL }, { type: 'graphene', Vg }, { n: Lo, d: dL }, ...mirror], 0, 's', gp);
      const sig = grapheneSigma(EV_NM / wl, mu_eV, T_eV, Gamma_eV);
      const em = solveTMM([wl], 1, 1.5,
        [...mirror, { n: Lo, d: dL }, { n: grapheneEffectiveIndex(sig, wl), d: 0.34 }, { n: Lo, d: dL }, ...mirror], 0, 's');
      maxDiff = Math.max(maxDiff, Math.abs(bc.A[0] - em.A[0]));
    }
    out['4. max|A_BC − A_effmed| [<0.01]'] = maxDiff.toExponential(2);
  }
  { // 5. Pauli blocking at λ = 4 µm (ħω = 0.31 eV)
    const mk = (Vg) => [{ n: air, d: 100 }, { type: 'graphene', Vg }, { n: air, d: 100 }];
    out['5a. A(4µm), Vg=8V unblocked [≈0.021]'] = solveTMM([4000], 1, 1, mk(8), 0, 's').A[0].toFixed(4);
    out['5b. A(4µm), Vg=100V blocked [≈0.002]'] = solveTMM([4000], 1, 1, mk(100), 0, 's').A[0].toFixed(4);
  }
  { // 6. no negative absorption at 45°, both pols, graphene inside the Bragg stack
    const layers = defaultLayers().map((l) => ({ n: { re: l.n, im: l.k }, d: l.d }));
    layers.splice(4, 0, { type: 'graphene', Vg: 30 });
    const wl = linspace(300, 900, 200);
    let minA = 1;
    for (const pol of ['s', 'p']) {
      const r = solveTMM(wl, 1, 1.5, layers, 45, pol);
      for (const a of r.A) minA = Math.min(minA, a);
    }
    out['6. min A @45° [≥ 0]'] = minA.toExponential(2);
  }

  console.table(out);
  return out;
}

// MARK: Control wiring
// ------------------------------------------------------------------------------------------

document.getElementById('calc-btn').addEventListener('click', calculate);
document.getElementById('add-layer-btn').addEventListener('click', addLayer);
document.getElementById('reset-btn').addEventListener('click', resetToExample);

document.getElementById('probe-slider').addEventListener('input', (e) => {
  document.getElementById('probe-value').textContent = `${Math.round(e.target.value)} nm`;
  drawChart();
  drawDiagram();
  updateReadout();
});

document.getElementById('vg-run-btn').addEventListener('click', runLayerVgSweep);
document.getElementById('vg-probe-slider').addEventListener('input', (e) => {
  const idx = parseInt(e.target.value);
  drawVgChart();
  updateVgReadout(idx);
});
drawVgChart();

renderLayers();
calculate();