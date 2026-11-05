/* ===== Synthetic ECG engine =====
   Coordinates: x = patient's left, y = inferior, z = anterior.
   Each wave component is a dipole with a 3D direction; lead voltage = amplitude * (dir · leadVector) * shape(t).
*/
const ECG = (() => {
  const n = v => { const m = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / m, v[1] / m, v[2] / m]; };
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const rad = d => d * Math.PI / 180;
  const LEADS = {
    I: [1, 0, 0], II: n([0.5, 0.866, 0]), III: n([-0.5, 0.866, 0]),
    aVR: n([-0.866, -0.5, 0]), aVL: n([0.866, -0.5, 0]), aVF: [0, 1, 0],
    V1: n([-0.45, 0.15, 0.88]), V2: n([-0.1, 0.18, 0.98]), V3: n([0.3, 0.2, 0.93]),
    V4: n([0.62, 0.22, 0.75]), V5: n([0.86, 0.2, 0.45]), V6: n([0.97, 0.18, 0.12])
  };
  const GAIN = { I: 1, II: 1, III: 1, aVR: 1, aVL: 1, aVF: 1, V1: 1.05, V2: 1.8, V3: 2.0, V4: 2.0, V5: 1.7, V6: 1.3 };

  function rng(seed) { let s = (seed || 7) >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  // rotate the frontal-plane (x,y) part of a vector by deg
  function rotF(v, deg) { const a = rad(deg), c = Math.cos(a), s = Math.sin(a); return [v[0] * c - v[1] * s, v[0] * s + v[1] * c, v[2]]; }
  const fdir = (deg, z = 0) => [Math.cos(rad(deg)), Math.sin(rad(deg)), z];

  /* ---- morphology templates (times in ms from QRS onset) ---- */
  function qrsAtoms(m) {
    const ax = m.axis ?? 55, d = ax - 35; // template's native axis ≈ 35°
    const A = (t, s, amp, dir, sR) => ({ t, s, sR: sR || s, amp, dir: n(dir) });
    let atoms;
    switch (m.type) {
      case 'rbbb':
        atoms = [A(14, 8, .28, [-.5, .3, .8]), A(36, 10, .95, rotF([.35, .75, .5], d)), A(52, 10, 1.05, rotF([.85, .3, -.45], d)),
          A(112, 19, 1.05, [-.75, .05, .75])];
        break;
      case 'lbbb':
        atoms = [A(22, 12, .45, [.75, .15, -.35]), A(62, 20, 1.45, rotF([.75, .35, -.5], d)), A(112, 17, 1.0, rotF([.85, .12, -.5], d))];
        break;
      case 'lafb': // initial inferior-rightward then superior-leftward
        atoms = [A(14, 8, .35, [-.15, .9, .35]), A(42, 11, 1.25, fdir(m.axis ?? -60, -.2)), A(58, 10, .7, fdir((m.axis ?? -60) + 15, -.45)), A(50, 11, .55, [.9, .3, -.25]), A(78, 8, .2, [-.3, -.4, -.6])];
        break;
      case 'lpfb':
        atoms = [A(14, 8, .35, [.6, -.7, .2]), A(42, 11, 1.25, fdir(m.axis ?? 115, .15)), A(58, 10, .7, fdir((m.axis ?? 115) - 10, -.45)), A(78, 8, .2, [-.3, -.4, -.6])];
        break;
      case 'pvc': // wide ectopic ventricular beat
        atoms = [A(45, 26, 1.3, m.dir1 || [-.55, -.75, .35], 30), A(110, 26, 1.0, m.dir2 || [-.7, -.55, .45], 30)];
        break;
      case 'vesc': // ventricular escape (wide, LBBB-like from RV)
        atoms = [A(50, 24, 1.1, [.7, .45, -.6], 28), A(115, 24, .75, [.8, .2, -.55], 28)];
        break;
      default:
        atoms = [A(14, 8.5, .27, [-.5, .3, .8]), A(37, 10.5, .95, rotF([.35, .75, .22], d)), A(53, 10.5, 1.1, rotF([.85, .3, -.62], d)), A(74, 8.5, .32, rotF([-.2, -.55, -.7], d))];
    }
    if (m.q) atoms.push(A(16, m.q.w || 13, m.q.amp, m.q.dir)); // infarct Q vector (points away from dead wall)
    if (m.lvh) atoms.forEach(a => { if (a.amp > .6) a.amp *= m.lvh; });
    if (m.widen) atoms.forEach(a => { a.t *= m.widen; a.s *= m.widen; a.sR *= m.widen; });
    if (m.scale) atoms.forEach(a => a.amp *= m.scale);
    return atoms;
  }
  function qrsEnd(m) { return ({ rbbb: 140, lbbb: 150, pvc: 160, vesc: 165 }[m.type] || 92) * (m.widen || 1); }
  function tAtoms(m, qt) {
    let dir, amp = .34, sL = 55, sR = 34;
    switch (m.type) {
      case 'rbbb': dir = [.6, .65, -.45]; break;            // inverted V1–V2
      case 'lbbb': dir = [-.65, -.2, .65]; amp = .38; break; // discordant
      case 'pvc': dir = m.tdir || [.6, .7, -.4]; amp = .6; sL = 50; sR = 40; break;
      case 'vesc': dir = [-.65, -.3, .6]; amp = .45; sL = 55; sR = 42; break;
      default: dir = [.55, .7, .45];
    }
    if (m.tdir) dir = m.tdir;
    if (m.tamp) amp = m.tamp;
    if (m.peakedT) { amp = 1.05; sL = 26; sR = 22; }
    const tp = qt - 2.4 * sR;
    return [{ t: tp, s: sL, sR, amp: amp * (m.scale && m.type !== 'pvc' ? Math.max(m.scale, .5) : 1), dir: n(dir) }];
  }
  function pAtoms(axis, kind, sc) {
    sc = sc ?? 1;
    if (kind === 'retro') return [{ t: 40, s: 18, sR: 18, amp: .2 * sc, dir: n([-.1, -1, .15]) }];
    const ra = n(rotF([.25, .9, .45], axis - 65)), la = n(rotF([.85, .45, -.4], axis - 65));
    return [{ t: 35, s: 17, sR: 17, amp: .14 * sc, dir: ra }, { t: 62, s: 18, sR: 18, amp: .13 * sc, dir: la }];
  }

  /* ---- rhythm scheduling: returns {atoms, extras} with absolute times (ms) ---- */
  function build(spec) {
    const R = rng(spec.seed || 11), atoms = [], dur = spec.dur || 10000;
    const morph = Object.assign({ type: 'normal' }, spec.qrs || {});
    const pr = spec.pr ?? 160, rate = spec.rate || 75, rr = 60000 / rate;
    const qt = spec.qt || Math.round(400 * Math.sqrt(rr / 1000) + 10);
    const pAx = spec.pAxis ?? 60;
    const add = (base, list, extra) => list.forEach(a => atoms.push(Object.assign({}, a, { t: a.t + base }, extra || {})));
    let qrsCount = 0;
    const addQRS = (t0, m, k) => {
      const mm = Object.assign({}, m);
      if (spec.alternans && (k % 2)) mm.scale = (mm.scale || 1) * spec.alternans;
      add(t0, qrsAtoms(mm));
      add(t0, tAtoms(mm, mm.type === 'pvc' ? qt + 40 : qt));
      if (spec.st && mm.type !== 'pvc') atoms.push({ plateau: true, t0: t0 + qrsEnd(mm) - 30, t1: t0 + qt - 90, amp: spec.st.mm / 10, dir: n(spec.st.dir) });
      qrsCount++;
    };
    const addP = (t0, ax, kind) => add(t0, pAtoms(ax, kind, spec.pScale));
    const R0 = spec.start ?? 260;
    const rh = spec.rhythm || 'sinus';

    if (rh === 'sinus' || rh === 'ectopicAtrial') {
      const ect = spec.ectopy || [];
      let pT = R0 - pr * 0.4, k = 0, cycle = 0;
      const cond = spec.conduction || 'normal';
      const prSeq = spec.prSeq || [pr];
      let wk = 0;
      let skipNext = false;
      while (pT < dur) {
        const ectHere = ect.find(e => e.beat === k);
        // sinus arrhythmia jitter
        const jitter = spec.jitter ? (R() - .5) * spec.jitter : 0;
        addP(pT, pAx, 'normal');
        let conducted = true, thisPR = pr;
        if (cond === 'mobitz1') { if (wk >= prSeq.length) { conducted = false; wk = 0; } else { thisPR = prSeq[wk]; wk++; } }
        else if (cond === 'mobitz2') { const nn = spec.ratio || 3; if ((k % nn) === nn - 1) conducted = false; }
        else if (cond === '2to1') { if (k % 2 === 1) conducted = false; }
        if (skipNext) { conducted = false; skipNext = false; }
        if (conducted) addQRS(pT + thisPR, morph, k);
        let next = pT + rr + jitter;
        if (ectHere) {
          const qrsT = pT + thisPR;
          if (ectHere.type === 'pvc') {
            const pv = Object.assign({ type: 'pvc' }, ectHere.morph || {});
            addQRS(qrsT + ectHere.coupling, pv, k);
            skipNext = true; // next sinus P is buried & not conducted → compensatory pause
          } else if (ectHere.type === 'pac') {
            const pacP = qrsT + ectHere.coupling - pr;
            addP(pacP, ectHere.pAxis ?? 10, 'normal');
            addQRS(pacP + (ectHere.pr || pr), morph, k);
            next = pacP + rr; // resets SA node → non-compensatory
          }
        }
        pT = next; k++;
      }
    } else if (rh === 'junctional') {
      for (let t = R0; t < dur; t += rr) { addQRS(t, morph, qrsCount); addP(t + (spec.retroDelay ?? 50), 0, 'retro'); }
    } else if (rh === 'chb') {
      const prr = 60000 / (spec.atrialRate || 88);
      for (let t = spec.pStart ?? 120; t < dur; t += prr) addP(t, pAx, 'normal');
      const esc = Object.assign({ type: spec.escape === 'junctional' ? 'normal' : 'vesc' }, spec.escapeMorph || {});
      for (let t = spec.vStart ?? 520; t < dur; t += rr) addQRS(t, esc, qrsCount);
    } else if (rh === 'af') {
      let times = [];
      if (spec.count) {
        const w = []; for (let i = 0; i < spec.count; i++) w.push(.55 + R() * .9);
        const span = dur - 520; const sw = w.slice(1).reduce((a, b) => a + b, 0);
        let t = 220; times.push(t); for (let i = 1; i < spec.count; i++) { t += w[i] / sw * span; times.push(t); }
      } else { let t = 200; while (t < dur) { times.push(t); t += rr * (.55 + R() * .9); } }
      times.forEach(t => addQRS(t, morph, qrsCount));
    } else if (rh === 'flutter') {
      const fr = 60000 / (spec.atrialRate || 300), ratio = spec.ratio || 2;
      let i = 0; for (let t = 60; t < dur; t += fr, i++) if (i % ratio === 1) addQRS(t + 120, morph, qrsCount);
    } else if (rh === 'vt') {
      for (let t = R0; t < dur; t += rr) addQRS(t, Object.assign({ type: 'pvc' }, spec.qrs || {}), qrsCount);
    }
    return { atoms, spec, dur, R: rng((spec.seed || 11) + 99) };
  }

  /* ---- sampling ---- */
  function sample(model, lead, t0, t1, step) {
    const L = LEADS[lead], g = GAIN[lead], sp = model.spec, out = [];
    const atoms = model.atoms.map(a => Object.assign({}, a, { p: dot(a.dir, L) }));
    const fv = sp.rhythm === 'af' ? n(sp.fDir || [.2, .6, .75]) : sp.rhythm === 'flutter' ? n([.1, -1, .35]) : null;
    const fP = fv ? dot(fv, L) : 0;
    const lowV = sp.lowVoltage || 1;
    for (let t = t0; t <= t1; t += step) {
      let v = 0;
      for (const a of atoms) {
        if (a.plateau) {
          if (t > a.t0 - 40 && t < a.t1 + 120) {
            const up = 1 / (1 + Math.exp(-(t - a.t0) / 9)), dn = 1 / (1 + Math.exp((t - a.t1) / 30));
            v += a.amp * a.p * up * dn;
          }
          continue;
        }
        const dt = t - a.t; if (dt < -4 * a.s || dt > 4 * a.sR) continue;
        const s = dt < 0 ? a.s : a.sR;
        v += a.amp * a.p * Math.exp(-dt * dt / (2 * s * s));
      }
      v *= lowV;
      if (sp.rhythm === 'af') {
        const am = sp.fAmp ?? .06;
        v += am * fP * (Math.sin(t * .0336 + 1.3) * .5 + Math.sin(t * .0431 + .2) * .35 + Math.sin(t * .0517 + 2.1) * .3) * (0.75 + 0.25 * Math.sin(t * .0011));
        v += am * .25 * Math.sin(t * .061 + (lead.length));
      }
      if (sp.rhythm === 'flutter') {
        const fr = 60000 / (sp.atrialRate || 300), ph = ((t - 60) % fr + fr) % fr / fr;
        const saw = ph < .72 ? (ph / .72) : (1 - (ph - .72) / .28); // slow rise toward peak of negative wave then quick return
        v += -(sp.fAmp ?? .22) * (saw - .5) * fP * -1;
      }
      v *= g;
      v += .018 * Math.sin(t * .0016 + lead.charCodeAt(0)) + (sp.noise ?? .006) * (model.R() - .5);
      out.push(v);
    }
    return out;
  }

  /* ---- SVG rendering on ECG paper ---- */
  const PX = 4; // px per mm
  function grid(w, h) {
    let s = '';
    for (let x = 0; x <= w; x += PX) s += `<line x1="${x}" y1="0" x2="${x}" y2="${h}" stroke="${(x / PX) % 5 ? '#f6cfcf' : '#e79a9a'}" stroke-width="${(x / PX) % 5 ? .5 : 1}"/>`;
    for (let y = 0; y <= h; y += PX) s += `<line x1="0" y1="${y}" x2="${w}" y2="${y}" stroke="${(y / PX) % 5 ? '#f6cfcf' : '#e79a9a'}" stroke-width="${(y / PX) % 5 ? .5 : 1}"/>`;
    return s;
  }
  function path(vals, x0, y0, step) {
    const mmPerMs = 25 / 1000, mmPerMv = 10;
    let d = '';
    vals.forEach((v, i) => { const x = x0 + i * step * mmPerMs * PX, y = y0 - v * mmPerMv * PX; d += (i ? 'L' : 'M') + x.toFixed(1) + ' ' + y.toFixed(1); });
    return `<path d="${d}" fill="none" stroke="#111" stroke-width="1.15" stroke-linejoin="round"/>`;
  }
  const cal = (x, y) => `<path d="M${x} ${y} h${PX} v${-10 * PX} h${5 * PX} v${10 * PX} h${PX}" fill="none" stroke="#111" stroke-width="1.15"/>`;
  const lab = (x, y, t) => `<text x="${x}" y="${y}" font-size="13" font-weight="700" font-family="Helvetica,Arial,sans-serif" fill="#222">${t}</text>`;

  function render(spec) {
    const model = build(spec), step = 4;
    const layout = spec.layout || '12';
    const calW = 8 * PX;
    if (layout === '12') {
      const rows = [['I', 'aVR', 'V1', 'V4'], ['II', 'aVL', 'V2', 'V5'], ['III', 'aVF', 'V3', 'V6']];
      const rowH = 34 * PX, W = calW + 250 * PX, H = rowH * 4 + 6 * PX;
      let g = grid(W, H), tr = '';
      rows.forEach((r, ri) => {
        const y0 = rowH * ri + rowH * .55;
        tr += cal(0, y0);
        r.forEach((ld, ci) => {
          const t0 = ci * 2500, vals = sample(model, ld, t0, t0 + 2500 - step, step);
          tr += path(vals, calW + ci * 62.5 * PX, y0, step) + lab(calW + ci * 62.5 * PX + 6, rowH * ri + 16, ld);
          if (ci) tr += `<line x1="${calW + ci * 62.5 * PX}" y1="${y0 - 6 * PX}" x2="${calW + ci * 62.5 * PX}" y2="${y0 + 3 * PX}" stroke="#111" stroke-width="1"/>`;
        });
      });
      const rl = spec.rhythmLead || 'II', y0 = rowH * 3 + rowH * .55;
      tr += cal(0, y0) + path(sample(model, rl, 0, 10000 - step, step), calW, y0, step) + lab(calW + 6, rowH * 3 + 16, rl);
      return svgWrap(W, H, g + tr);
    } else { // rhythm strip(s)
      const leads = spec.leads || ['II'], rowH = (spec.rowH || 30) * PX, secs = (spec.dur || 10000);
      const W = calW + secs / 1000 * 25 * PX, H = rowH * leads.length + 2 * PX;
      let g = grid(W, H), tr = '';
      leads.forEach((ld, i) => {
        const y0 = rowH * i + rowH * .6;
        tr += cal(0, y0) + path(sample(model, ld, 0, secs - step, step), calW, y0, step) + lab(calW + 6, rowH * i + 16, ld);
      });
      return svgWrap(W, H, g + tr);
    }
  }
  function svgWrap(W, H, inner) {
    return `<svg class="ecgsvg" viewBox="0 0 ${W} ${H + 18}" width="${W}" height="${H + 18}" xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMinYMin meet"><rect width="${W}" height="${H + 18}" fill="#fff8f8"/>${inner}<text x="6" y="${H + 13}" font-size="11" font-family="Helvetica,Arial,sans-serif" fill="#7a4a4a">25 mm/s · 10 mm/mV · small box = 40 ms / 0.1 mV · large box = 200 ms / 0.5 mV</text></svg>`;
  }
  return { render, build, sample, LEADS };
})();
