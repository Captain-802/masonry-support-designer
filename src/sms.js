/* Calculation report for the steel masonry support (BS 5950).
   Builds the calculation sheet line for line in the reference layout, with its
   "Steel beam torsion design" v2.0.02 part: same headings, labels, formulas, number formats and check
   messages (typos and spacing included), taken from 26 the reference outputs (20 GUI runs, 6 worked examples).

   build(p, o, opts) -> items       p, o: SMS.calc input and result (units N, mm)
     item.t : 'title' | 'sub' | 'ver' | 'sketch' | 'h' (bold heading) | 'row' | 'row2' (expression at the
              second tab stop) | 'txt' (plain paragraph) | 'chk' (bold italic, right aligned)
     item.l : label (markup)          item.e : expression (markup)          item.s : text (markup)
   Markup: _{..} subscript, ^{..} superscript, **..** bold (values); Greek letters, × and √ are drawn in a
   symbol font. plain(markup) gives the text exactly as the reference copies it to the clipboard (Symbol-font
   letters, sub/superscripts flattened), which the tests compare with the reference outputs.
   Lines marked UNVERIFIED have no the reference sample yet (see tests). */
(function (root) {
  'use strict';
  const PI_T = 3.1415926536;   // the reference value of pi (see sms.js)
  // ------------------------------------------------------------------ text helpers
  const SYMBOL = { 'α': 'a', 'β': 'b', 'γ': 'g', 'δ': 'd', 'ε': 'e', 'η': 'h', 'θ': 'q', 'λ': 'l', 'μ': 'm',
    'π': 'p', 'ρ': 'r', 'σ': 's', 'τ': 't', 'φ': 'f', 'χ': 'c', 'ψ': 'y', 'ω': 'w', 'Δ': 'D',
    '×': '´', '√': 'Ö', '≤': '£', '≥': '³' };
  const plain = m => String(m).replace(/\*\*/g, '').replace(/[_^]\{([^}]*)\}/g, '$1').replace(/[^\x00-\x7f]/g, ch => SYMBOL[ch] || ch);

  // Rounding of printed values: scale by 10^d in double arithmetic, then round half away from zero.
  // So 2550 N/m prints as 2.6 kN/m (2.55 x 10 = 25.5) while 0.012499999999999997 m
  // prints as 12 mm; JavaScript's toFixed would give 2.5 and 12. "-0" prints as "0".
  function tround(x, d) {
    const k = 10 ** d, v = x * k, i = Math.trunc(v), fr = v - i;
    return (Math.abs(fr) >= 0.5 ? i + Math.sign(fr) : i) / k;
  }
  function fx(x, d) {
    if (!Number.isFinite(x)) return String(x);
    let s = tround(x, d).toFixed(d);
    if (/^-0(\.0*)?$/.test(s)) s = s.slice(1);
    return s;
  }
  // 3 significant figures, scientific (2.64×10^-2) or engineering (14.0×10^-3), bold as the reference prints them
  function sci(x) {
    if (x === 0 || !Number.isFinite(x)) return `**${fx(x, 2)}**`;
    let e = Math.floor(Math.log10(Math.abs(x))), m = x / 10 ** e;
    if (Math.abs(tround(m, 2)) >= 10) { e += 1; m = x / 10 ** e; }
    return `**${fx(m, 2)}**×**10^{${e}}**`;
  }
  function eng(x) {
    if (x === 0 || !Number.isFinite(x)) return `**${fx(x, 2)}**`;
    let e = 3 * Math.floor(Math.floor(Math.log10(Math.abs(x))) / 3), m = x / 10 ** e;
    let d = Math.abs(m) < 10 ? 2 : Math.abs(m) < 100 ? 1 : 0;
    if (Math.abs(tround(m, d)) >= 1000) { e += 3; m = x / 10 ** e; d = 2; }
    return `**${fx(m, d)}${d === 0 ? '.' : ''}**×**10^{${e}}**`;
  }
  const B = (x, d) => `**${fx(x, d)}**`;

  // Table 13 wording: both_free (worked examples), both_full (live run V1), comp_full (V2); the partial rows follow the table
  const RESTRAINT = {
    both_free: ['Compression flanges laterally restrained', 'Both flanges free to rotate on plan'],
    both_full: ['Compression flanges laterally restrained', 'Both flanges fully restrained against rotation on plan'],
    comp_full: ['Compression flanges laterally restrained', 'Compression flange fully restrained against rotation on plan'],
    both_partial: ['Compression flanges laterally restrained', 'Both flanges partially restrained against rotation on plan'],
    comp_partial: ['Compression flanges laterally restrained', 'Compression flange partially restrained against rotation on plan'],
  };
  const K_TABLE13 = { both_full: [0.7, 0.85], comp_full: [0.75, 0.9], both_partial: [0.8, 0.95], comp_partial: [0.85, 1.0], both_free: [1.0, 1.2] };
  function restraintOf(p) {
    if (p.restraint && RESTRAINT[p.restraint]) return p.restraint;
    for (const col of p.destab ? [1, 0] : [0, 1]) for (const [k, v] of Object.entries(K_TABLE13)) if (Math.abs(v[col] - p.K_LT) < 1e-9) return k;
    return 'both_free';
  }
  const kText = k => (Math.round(k * 10) === k * 10 ? k.toFixed(1) : String(+k.toFixed(3)));

  // ------------------------------------------------------------------ SI shadow
  // the reference evaluates each printed formula in SI units (m, N, Pa) in the order printed and rounds that double.
  // The short formulas are re-evaluated here the same way, so values on a rounding tie print as the reference prints
  // them (e.g. 0.1025 m -> "102" mm, 1112649.9999999998 N -> "1112.6" kN). Long chains use the engine values.
  function si(p, o) {
    const S = {}, mb = p.mb || null, tor = p.design_method === 'Masonry support and torsion', m = p.sb_method;
    const bolt = p.fixing === 'Bolt', slabS = p.sb_load === 'Slab' && m !== 'SCI', slabM = tor && p.mb_load === 'Slab';
    const E = 205e9, RHO_S = 78500, gG = p.gG, gQ = p.gQ, mm = x => x / 1000;
    const c = mm(p.c), emb = tor ? mm(p.e_mb || 0) : 0, s = bolt ? mm(p.shim || 0) : 0, e1 = mm(p.e1 || 0);
    const bmsb = mm(p.b_msb), hmsb = mm(p.h_msb), rho = p.rho_m_sb * 1e9;
    const PGa = p.PGadd_sb * 1000, PQa = p.PQadd_sb * 1000;
    const pysb = o.pysb * 1e6;
    let lh, lv, tsb;
    if (m === 'Plate') { lh = mm(p.lh); tsb = mm(p.tsb); }
    else if (m === 'Cantilever') { const long = (p.orientation || 'Long') === 'Long'; lh = mm(long ? p.angle.D : p.angle.B); lv = mm(long ? p.angle.B : p.angle.D); tsb = mm(p.angle.t); }
    else { lh = mm(p.lh); lv = mm(p.lv); tsb = mm(p.tsb); }
    Object.assign(S, { lh, lv, tsb, c, emb, s, e1 });
    // loads on the support member
    S.P1 = slabS ? p.PG_ssb * 1000 * gG + p.PQ_ssb * 1000 * gQ : (bmsb * hmsb * rho + PGa) * gG + PQa * gQ;
    S.P1SLS = slabS ? p.PG_ssb * 1000 + p.PQ_ssb * 1000 : bmsb * hmsb * rho + PGa + PQa;
    S.dm = m === 'Plate' ? lh + emb - c : lh + s + emb - c;
    const Zp = 1 * tsb ** 2 / (6 * 1);
    if (m === 'Plate') {
      const lplate = mm(p.lplate), Bmb = mm(mb.B), Dmb = mm(mb.D), Amb = mb.A / 1e6, Ixmb = mb.Ix / 1e12;
      S.Asbu = tsb * lplate; S.cyysb = lh / 2 - (lplate - lh) / 2;
      S.ye = (Dmb + tsb) * S.Asbu / (2 * (Amb + S.Asbu));
      S.Ixx = (Ixmb + Amb * S.ye ** 2) + S.Asbu * (Dmb / 2 + tsb / 2 - S.ye) ** 2;
      S.Zxx = S.Ixx / (Dmb / 2 + tsb - S.ye);
      S.Zp = Zp; S.Mxp = S.P1 * e1; S.Mc = 1.2 * Zp * pysb; S.als = S.Mxp / S.Mc;
      S.am = e1; S.bm = lh - e1; S.al = S.am / (S.am + S.bm); S.Ieff = tsb ** 3 / 12;
      S.delta = (S.al ** 2 * (3 - S.al) / 6) * (S.P1SLS * (S.am + S.bm) ** 3) / (E * S.Ieff);
      if (p.con_stage) {
        S.e1c = c - emb + bmsb / 2 - s; S.P1c = 1.2 * bmsb * rho * gG; S.Mxpc = S.P1c * S.e1c; S.alsc = S.Mxpc / S.Mc;
        S.P1cS = 1.2 * bmsb * rho; S.amc = S.e1c; S.bmc = lh - S.e1c; S.alc = S.amc / (S.amc + S.bmc);
        S.deltac = (S.alc ** 2 * (3 - S.alc) / 6) * (S.P1cS * (S.amc + S.bmc) ** 3) / (E * S.Ieff);
      }
      const sweld = mm(p.s_weld), L = mm(p.L);
      S.aw = 1 / Math.sqrt(2) * sweld;
      S.RA = S.P1 * Math.max((1 + (3 * e1) / (2 * Bmb / 2)), 1.4);
      S.Rp = (lh + Bmb) * tsb * pysb; S.Rl = 2 * S.Rp / L; S.Rh = S.P1 * e1 / (sweld / 2 + tsb / 2);
      S.Rweld = (S.RA ** 2 + S.Rl ** 2 + S.Rh ** 2) ** 0.5; S.pweld = o.pweld * 1e6; S.pcw = S.aw * S.pweld;
      S.dlim = c > 0 ? Math.min((1 + S.dm / c) * 0.001, 0.002) : 0.002;
    } else if (m === 'Cantilever') {
      S.Asb = p.angle.A / 1e6; S.cyysb = p.angle.cyy / 1000;
      const lb = bolt ? lv - mm(p.eb) : lv, hb = (o.hb || 0) / 1000, Bb = bolt ? mm(p.Bbolt) : 0;
      S.lb = lb; S.l = lb + hb; S.M1 = S.P1 * e1; S.Zp = Zp; S.Mc = 1.2 * pysb * Zp;
      if (bolt) {
        S.lamm = Bb / (lh + lb); S.Bstr = Math.min(2 * Bb / (0.4 + S.lamm), Bb); S.sigma = Bb * S.M1 / (S.Bstr * Zp);
        S.Bdef = Math.min(2 * Bb / (0.8 + S.lamm), Bb); S.Ieff = S.Bdef * tsb ** 3 / 12;
      } else { S.sigma = S.M1 / Zp; S.Ieff = 1 * tsb ** 3 / (12 * 1); }
      // a_m, b_m and d_m of angles come out of the reference as exact decimals (worked examples X01-X04): mm arithmetic
      const long = (p.orientation || 'Long') === 'Long', lhm = long ? p.angle.D : p.angle.B;
      S.am = (p.e1 + (o.hb || 0)) / 1000; S.bm = (lhm - p.e1) / 1000; S.dm = (lhm + (bolt ? p.shim || 0 : 0) + (tor ? p.e_mb || 0 : 0) - p.c) / 1000;
      S.pdef = S.am / S.l; S.qdef = S.bm / S.l;
      const pq = 4 * S.pdef ** 2 + 6 * S.pdef * S.qdef + 3 * S.pdef + 3 * S.qdef;
      S.delta = bolt ? (Bb * S.P1SLS * S.l ** 3 * S.pdef) / (12 * E * S.Ieff) * pq : (S.P1SLS * S.l ** 3 * S.pdef) / (12 * E * S.Ieff) * pq;
      S.dlim = Math.min((1 + S.dm / c) * 0.001, 0.002);
      if (!bolt) { S.aw = 1 / Math.sqrt(2) * mm(p.s_weld); S.fweld = Math.sqrt(S.P1 ** 2 + (S.M1 / S.l) ** 2); S.Pweld = S.fweld / (p.l_weld * S.aw); }
      if (bolt) { S.Fs = S.P1 * Bb; S.Ft = Bb * S.M1 / S.l; }
    } else {
      const r1 = mm(p.r1sb), Bb = mm(p.Bbolt), lb = lv - mm(p.eb);
      S.Asbu = tsb * (lh + lv + PI_T * (tsb / 2 + r1) / 2);
      S.cyysb = (lv * tsb ** 2 / 2 + lh * tsb * (lh - tsb) / 2) / ((lh + lv) * tsb - tsb ** 2);
      S.lb = lb; S.l = lb - (tsb + r1); S.g = lv - lb;
      S.e = Math.max(0.01, Math.min(bmsb / 2, Math.min(((16.5 * tsb / (S.l + S.g)) - 0.73), 1) * S.dm));
      S.am = c - emb - (s + tsb + r1) + PI_T * (tsb / 2 + r1) / 2 + S.e; S.bm = S.dm - S.e;
      S.lamm = Bb / (S.l + S.am + S.bm); S.Bstr = Math.min(2 * Bb / (0.4 + S.lamm), Bb);
      S.sigma = 6 * S.P1 * Bb * (c - emb + S.e - s - 0.5 * tsb) / (S.Bstr * tsb ** 2);
      S.Bdef = Math.min(2 * Bb / (0.8 + S.lamm), Bb); S.Ieff = S.Bdef * tsb ** 3 / 12; S.pdef = S.am / S.l; S.qdef = S.bm / S.l;
      S.delta = (Bb * S.P1SLS * S.l ** 3 * S.pdef) / (12 * E * S.Ieff) * (4 * S.pdef ** 2 + 6 * S.pdef * S.qdef + 3 * S.pdef + 3 * S.qdef);
      S.dlim = Math.min((1 + S.dm / c) * 0.001, 0.002);
      S.Fs = S.P1 * Bb; S.Ft = Bb * S.P1 * S.am / S.l;
    }
    if (bolt) S.Fcomb = S.Fs / o.Ps_bolt + S.Ft / o.Pt_bolt;
    if (!tor) return S;
    // torsion beam loads, eccentricities and torque
    const Bmb = mm(mb.B), tmb = mm(mb.t), L = mm(p.L), e0 = (o.e0mb || 0) / 1000;
    const Asw = m === 'Cantilever' ? S.Asb : S.Asbu;
    S.w1 = S.P1; S.w1S = S.P1SLS;
    if (slabM) { S.w2 = p.PG_smb * 1000 * gG + p.PQ_smb * 1000 * gQ; S.w2S = p.PG_smb * 1000 + p.PQ_smb * 1000; }
    else {
      const bmmb = mm(p.b_mmb), hmmb = mm(p.h_mmb), rmb = p.rho_m_mb * 1e9, PGm = p.PGadd_mb * 1000, PQm = p.PQadd_mb * 1000;
      S.w2 = (hmmb * bmmb * rmb + PGm) * gG + PQm * gQ; S.w2S = hmmb * bmmb * rmb + PGm + PQm;
    }
    S.w3 = Asw * RHO_S * gG; S.w3S = Asw * RHO_S;
    const isC = mb.type === 'C', face = isC ? tmb / 2 : Bmb / 2;
    S.e1mb = slabS ? (isC ? e1 + s + face - e0 : e1 - e0 + s + face) : isC ? c + (bmsb + tmb) / 2 - emb - e0 : (Bmb + bmsb) / 2 + c - emb;
    S.e2mb = slabM ? (isC ? -mm(p.e2) - e0 + face : Bmb / 2 - mm(p.e2)) : isC ? -(mm(p.b_mmb) - tmb) / 2 - emb - e0 : (Bmb - mm(p.b_mmb)) / 2 - emb;
    S.e3mb = isC ? S.cyysb - e0 + tmb / 2 : Bmb / 2 + S.cyysb;
    if (p.mirror && isC) { const sh = (o.mirrorShift || 0) / 1000; S.e1mb += sh; S.e2mb += sh; S.e3mb += sh; }
    S.TqULS = Math.abs(S.w1 * S.e1mb + S.w2 * S.e2mb + S.w3 * S.e3mb); S.Tq = S.TqULS * L;
    S.TqSLS = Math.abs(S.w1S * S.e1mb + S.w2S * S.e2mb + S.w3S * S.e3mb); S.Tqu = S.TqSLS * L;
    // beam actions
    const Amb = mb.A / 1e6, py = o.py * 1e6;
    S.w4 = Amb * RHO_S * gG; S.w4S = Amb * RHO_S;
    S.wULS = S.w1 + S.w2 + S.w3 + S.w4; S.wSLS = S.w1S + S.w2S + S.w3S + S.w4S;
    S.Mx = S.wULS * L ** 2 / 8 + p.Mx_add / 1000; S.Fv = S.wULS * L / 2 + p.Fv_add;
    S.eps = Math.sqrt(275e6 / py); S.LE = mm(p.LLT) * p.K_LT;
    const Avy = mb.type === 'HR' ? (mb.A / 1e6) * mm(mb.D) / (mm(mb.D) + Bmb) : tmb * mm(mb.D);
    S.Pv = 0.6 * py * Avy;
    const plastic = o.Class === 'Plastic' || o.Class === 'Compact', Sx = mb.Sx / 1e9, Zx = mb.Zx / 1e9;
    S.Mcxu = py * (plastic ? Sx : Zx); S.Mcx = Math.min(S.Mcxu, 1.2 * py * Zx);
    if (o.highShear) {           // Cl. 4.2.5.3 (the reference labels it 4.2.5.1); S_vx as the engine
      const Svx = mb.type === 'HR' ? Avy * Avy / (8 * mm(mb.t)) : tmb * mm(mb.D) ** 2 / 4;
      S.rhoy = (2 * (S.Fv / S.Pv) - 1) ** 2; S.Mcxu = py * (Sx - S.rhoy * Svx); S.Mcx = Math.min(S.Mcxu, 1.2 * py * Zx);
    }
    if (mb.type !== 'HR') {
      const ry = mm(mb.ry);
      S.lam = S.LE / ry; S.v = 1 / (1 + 0.05 * (S.lam / mb.x) ** 2) ** 0.25; S.bw = plastic ? 1.0 : Zx / Sx;
      S.lamLT = mb.u * S.v * S.lam * Math.sqrt(S.bw); S.lamL0 = 0.4 * Math.sqrt(PI_T ** 2 * E / py);
      S.pE = PI_T ** 2 * E / S.lamLT ** 2; S.eta = Math.max(7.0 * (S.lamLT - S.lamL0) / 1000, 0);
      S.phiLT = (py + (S.eta + 1) * S.pE) / 2; S.pb = S.pE * py / (S.phiLT + Math.sqrt(S.phiLT ** 2 - S.pE * py));
      S.Mb = S.pb * (plastic ? Sx : Zx);
    } else {
      S.lam = S.LE / mm(mb.ry); S.Mb = S.Mcx;
      if (o.phib !== undefined) {   // beyond the Table 15 limit: Annex B.2.6.1 (live reference run V6)
        const Iy = mb.Iy / 1e12, Ix = mb.Ix / 1e12, J = mb.J / 1e12, A = mb.A / 1e6;
        S.gb = (1 - Iy / Ix) * (1 - J / (2.6 * Ix)); S.phib = Math.sqrt(Sx ** 2 * S.gb / (A * J));
        S.bw = plastic ? 1.0 : Zx / Sx; S.lamLT = 2.25 * Math.sqrt(S.phib * S.lam * S.bw);
        S.lamL0 = 0.4 * Math.sqrt(PI_T ** 2 * E / py); S.pE = PI_T ** 2 * E / S.lamLT ** 2;
        S.eta = Math.max(7.0 * (S.lamLT - S.lamL0) / 1000, 0); S.phiLT = (py + (S.eta + 1) * S.pE) / 2;
        S.pb = S.pE * py / (S.phiLT + Math.sqrt(S.phiLT ** 2 - S.pE * py)); S.Mb = S.pb * (plastic ? Sx : Zx);
      }
    }
    const phi = o.phi_deg * PI_T / 180;
    S.Myt = S.Mx * phi; S.sbyt = S.Myt / (mb.Zy / 1e9); S.sbx = S.Mx / Zx;
    S.dy = 5 * S.wSLS * L ** 4 / (384 * E * mb.Ix / 1e12) + p.d_add / 1000;
    S.dlimb = Math.min(L / p.k_delta, p.d_lim_abs / 1000);
    return S;
  }

  // ------------------------------------------------------------------ report
  function build(p, o, opts) {
    opts = opts || {};
    const items = [];
    const R = (l, e) => items.push({ t: 'row', l, e });
    const R2 = (l, e) => items.push({ t: 'row2', l, e });
    const H = s => items.push({ t: 'h', s });
    const T = s => items.push({ t: 'txt', s });
    const C = s => items.push({ t: 'chk', s });
    const ok = u => u <= 1 + 1e-12;
    const tor = p.design_method === 'Masonry support and torsion', m = p.sb_method, bolt = p.fixing === 'Bolt';
    const slabS = p.sb_load === 'Slab' && m !== 'SCI', slabM = tor && p.mb_load === 'Slab';
    // "Direct loads" (a page option, not a reference input): the typed dead and imposed line loads are listed in place of
    // the leaf height and density, then their totals G and Q carry on through the reference formulas
    const dirS = !slabS && !!p.direct_sb, dirM = tor && !slabM && !!p.direct_mb;
    const kNm = w => B(w, Math.abs(Math.round(w * 10) - w * 10) < 1e-9 ? 1 : 2);
    const short = s => (s.length > 40 ? s.slice(0, 39).trimEnd() + '…' : s);
    function loadList(lists, where, sym) {
      for (const [k, word] of [['G', 'dead'], ['Q', 'imposed']]) {
        const rows = (lists[k] || []).filter(r => Number.isFinite(r.w)), tot = rows.reduce((a, r) => a + r.w, 0);
        const lab = (r, i) => short(r.label || `Unfact ${word} load${rows.length > 1 ? ' ' + (i + 1) : ''} on ${where}`);
        if (rows.length === 1) { R(lab(rows[0], 0), `${k}_{${sym}} = ${kNm(tot)} kN/m`); continue; }
        rows.forEach((r, i) => R(lab(r, i), `${k}_{${sym},${i + 1}} = ${kNm(r.w)} kN/m`));
        R(`Total unfact ${word} load on ${where}`, `${k}_{${sym}} = ${rows.length ? rows.map((r, i) => `${k}_{${sym},${i + 1}}`).join(' + ') + ' = ' : ''}${kNm(tot)} kN/m`);
      }
    }
    const mb = p.mb, u = o.util_sb || {};
    const N2 = ' N/mm^{2}';
    const S = si(p, o), D = (x, f, d) => B(x / f, d), MM = x => B(x / 1000 / 0.001, 0);   // SI value / display unit

    items.push({ t: 'title', s: 'STEEL MASONRY SUPPORT' });
    items.push({ t: 'sub', s: 'In accordance with BS5950-1:2000 incorporating Corrigendum No.1' });
    items.push({ t: 'ver', s: opts.versionLine || '' });
    items.push({ t: 'sketch' });

    // ---- steel member details
    H('Steel member details');
    if (tor) R('Torsion beam', mb.name);
    let lh, lv;
    if (m === 'Plate') {
      lh = p.lh;
      R('Masonry support plate', 'User');
      R('Steel grade of support plate', p.grade_sb);
      R('Design strength of support plate', `p_{ysb} = ${B(o.pysb, 0)}${N2}`);
      R('Modulus of elasticity', 'E = 205000 N/mm^{2}');
      R('Constant', `ε = √(275N/mm^{2} / p_{ysb}) = ${D(Math.sqrt(275e6 / (o.pysb * 1e6)), 1, 3)}`);
      R('Length of plate beyond beam', `l_{h} = ${D(S.lh, 1e-3, 0)} mm`);
      R('Total length of plate', `l_{plate} = ${MM(p.lplate)} mm`);
      R('Thickness of plate', `t_{sb} = ${D(S.tsb, 1e-3, 0)} mm`);
      R('Width of main beam', `B_{mb} = ${MM(mb.B)} mm`);
      R('Area of plate', `A_{sbu} = t_{sb} × l_{plate} = ${D(S.Asbu, 1e-6, 1)} mm^{2}`);
      R('Distance from weld position to CoG', `c_{yysb} = l_{h} / 2 - (l_{plate} - l_{h}) / 2 = ${D(S.cyysb, 1e-3, 0)} mm`);
    } else if (m === 'Cantilever') {
      const long = (p.orientation || 'Long') === 'Long';
      lh = long ? p.angle.D : p.angle.B; lv = long ? p.angle.B : p.angle.D;
      R('Masonry support angle', p.angle.name);
      R('Steel grade of support angle', p.grade_sb);
      R('Design strength of support angle', `p_{ysb} = ${B(o.pysb, 0)}${N2}`);
      R('Modulus of elasticity', 'E = 205000 N/mm^{2}');
      R('Constant', `ε = √(275N/mm^{2} / p_{ysb}) = ${D(Math.sqrt(275e6 / (o.pysb * 1e6)), 1, 3)}`);
      R('Length of vertical leg', `l_{v} = ${long ? 'B' : 'D'}_{sb} = ${D(S.lv, 1e-3, 0)} mm`);
      R('Length of horizontal leg', `l_{h} = ${long ? 'D' : 'B'}_{sb} = ${D(S.lh, 1e-3, 0)} mm`);
      R('Orientation of angle', long ? 'Long leg horizontal' : 'Short leg horizontal');
      if (tor) R('Dist from soffit of angle to soffit of main beam', `h_{bf} = ${MM(p.h_bf || 0)} mm`);
    } else {
      lh = p.lh; lv = p.lv;
      R('Masonry support angle', 'User');
      R('Stainless steel grade of support angle', p.grade_sb);
      R('0.2% proof stress of support angle (SCI AD187)', `Y_{0.2} = ${B(o.pysb * (p.gmss || 1), 0)}${N2}`);
      R('Design strength of support angle', `p_{ysb} = Y_{0.2} / γ_{mss} = ${B(o.pysb, 0)}${N2}`);
      R('Length of vertical leg', `l_{v} = ${D(S.lv, 1e-3, 0)} mm`);
      R('Length of horizontal leg', `l_{h} = ${D(S.lh, 1e-3, 0)} mm`);
      R('Thickness of angle', `t_{sb} = ${D(S.tsb, 1e-3, 0)} mm`);
      R('Internal radius of angle', `r_{1sb} = ${D(p.r1sb / 1000, 1e-3, 1)} mm`);
      R('Area of angle', `A_{sbu} = t_{sb} × (l_{h} + l_{v} + π × (t_{sb} / 2 + r_{1sb}) / 2) = ${D(S.Asbu, 1e-6, 0)} mm^{2}`);
      R('Distance from back of vertical to CoG', `c_{yysb} =  (l_{v} × t_{sb}^{2} / 2 + l_{h} × t_{sb} × (l_{h} - t_{sb}) / 2) / ((l_{h} + l_{v}) × t_{sb} - t_{sb}^{2}) = ${D(S.cyysb, 1e-3, 0)} mm`);
      if (tor) R('Dist from soffit of angle to soffit of main beam', `h_{bf} = ${B(0, 0)} mm`);
    }

    // ---- supported materials
    H('Supported materials detail');
    if (tor) {
      if (slabM) {
        R('Depth of slab on main beam', `d_{smb} = ${MM(p.d_smb)} mm`);
        R('Unfact dead slab reaction on main beam', `P_{Gsmb} = ${B(p.PG_smb, 1)} kN/m`);
        R('Unfact imposed slab reaction on main beam', `P_{Qsmb} = ${B(p.PQ_smb, 1)} kN/m`);
        R('Eccentricity of main beam material', `e_{mb} = ${MM(p.e_mb)} mm`);
      } else if (dirM) {
        R('Width of masonry on main beam', `b_{mmb} = ${MM(p.b_mmb)} mm`);
        R('Eccentricity of main beam material', `e_{mb} = ${MM(p.e_mb)} mm`);
        loadList(p.direct_mb, 'main beam', 'mb');
      } else {
        R('Density of masonry on main beam', `ρ_{m,mb} = ${D(p.rho_m_mb * 1e9, 1000, 1)} kN/m^{3}`);
        R('Width of masonry on main beam', `b_{mmb} = ${MM(p.b_mmb)} mm`);
        R('Height of masonry on main beam', `h_{mmb} = ${MM(p.h_mmb)} mm`);
        R('Eccentricity of main beam material', `e_{mb} = ${MM(p.e_mb)} mm`);
        R('Add dead force main beam (not from masonry)', `P_{Gaddmb} = ${B(p.PGadd_mb, 1)} kN/m`);
        R('Add live force main beam (not from masonry)', `P_{Qaddmb} = ${B(p.PQadd_mb, 1)} kN/m`);
      }
    }
    if (slabS) {
      R('Depth of slab on support beam', `d_{ssb} = ${MM(p.d_ssb)} mm`);
      R('Unfact dead slab reaction on support beam', `P_{Gssb} = ${B(p.PG_ssb, 1)} kN/m`);
      R('Unfact imposed slab reaction on support beam', `P_{Qssb} = ${B(p.PQ_ssb, 1)} kN/m`);
    } else if (dirS) {
      R('Width of masonry on support beam', `b_{msb} = ${MM(p.b_msb)} mm`);
      loadList(p.direct_sb, 'support beam', 'sb');
    } else {
      R('Density of masonry on support beam', `ρ_{m,sb} = ${D(p.rho_m_sb * 1e9, 1000, 1)} kN/m^{3}`);
      R('Width of masonry on support beam', `b_{msb} = ${MM(p.b_msb)} mm`);
      R('Height of masonry on support beam', `h_{msb} = ${MM(p.h_msb)} mm`);
      if (m !== 'SCI') {
        R('Add dead force support beam (not from masonry)', `P_{Gaddsb} = ${B(p.PGadd_sb, 1)} kN/m`);
        R('Add live force support beam (not from masonry)', `P_{Qaddsb} = ${B(p.PQadd_sb, 1)} kN/m`);
      }
    }

    // ---- geometry
    H('Geometry');
    const hb = o.hb || 0;
    if (bolt) {
      R('Bolt spacing', `B_{bolt} = ${MM(p.Bbolt)} mm`);
      R('Distance from top of angle to bolt centreline', `e_{b} = ${MM(p.eb)} mm`);
      R('Dist bearing point of angle to bolt centreline', `l = l_{v} - e_{b}${hb < 0 ? ' + h_{bf}' : ''} = ${D(m === 'SCI' ? S.lb : S.l, 1e-3, 0)} mm`);
    } else if (m === 'Cantilever') R('Dist bearing point of angle to weld', `l = l_{v}${hb < 0 ? ' + h_{bf}' : ''} = ${D(S.l, 1e-3, 0)} mm`);   // h_bf form: live reference run V3
    R('Cavity width', `c = ${D(S.c, 1e-3, 0)} mm`);
    if (bolt) R('Shim thickness', `s = ${MM(p.shim || 0)} mm`);
    if (m === 'Plate') R('Supported width of masonry', `d_{m} = l_{h} + e_{mb} - c = ${D(S.dm, 1e-3, 0)} mm`);
    else R('Supported width of masonry', `d_{m} = l_{h} + s${tor ? ' + e_{mb}' : ''} - c = ${D(S.dm, 1e-3, 0)} mm`);

    const P1f = slabS ? 'P_{1} = P_{Gssb} × γ_{fG} + P_{Qssb} × γ_{fQ}' : dirS ? 'P_{1} = G_{sb} × γ_{fG} + Q_{sb} × γ_{fQ}' : 'P_{1} = (b_{msb} × h_{msb} × ρ_{m,sb} + P_{Gaddsb}) × γ_{fG} + P_{Qaddsb} × γ_{fQ}';
    const P1Sf = slabS ? 'P_{1SLS} = P_{Gssb} + P_{Qssb}' : dirS ? 'P_{1SLS} = G_{sb} + Q_{sb}' : 'P_{1SLS} = b_{msb} × h_{msb} × ρ_{m,sb} + P_{Gaddsb} + P_{Qaddsb}';
    const what = slabS ? 'slab' : dirS ? 'load' : 'masonry';

    if (m === 'Plate') {
      // ---- welded plate (SCI P110)
      H('Biaxial stress effects in the plate (SCI-P-110)');
      R('Maximum overall bending moment', `M_{x} = ${D(S.Mx, 1000, 1)} kNm`);
      R('Dist to NA combined section (CoG torsion beam)', `y_{e,all} = (D_{mb} + t_{sb}) × A_{sbu} / (2 × (A_{mb} + A_{sbu})) = ${D(S.ye, 1e-3, 0)} mm`);
      R('Second moment of area of combined section', `I_{xx,all} = (I_{xxmb} + A_{mb} × y_{e,all}^{2})+ A_{sbu} × (D_{mb} / 2 + t_{sb} / 2 - y_{e,all})^{2} = ${D(S.Ixx, 1e-8, 0)} cm^{4}`);
      R('Elastic section modulus of combined section', `Z_{xx,all} = I_{xx,all} / (D_{mb} / 2 + t_{sb} - y_{e,all}) = ${D(S.Zxx, 1e-6, 2)} cm^{3}`);
      R('Section modulus of plate', `Z_{xx,plate} = 1m × t_{sb}^{2} / (6 × 1m) = ${D(S.Zp, 1e-6, 2)} cm^{3}/m`);
      R(`Eccentricity of support beam ${what}`, `e_{1} = ${D(S.e1, 1e-3, 0)} mm`);
      R(`Force of ${what} on support plate`, `${P1f} = ${D(S.P1, 1000, 1)} kN/m`);
      R('Bending at heel', `M_{x,plate} = P_{1} × e_{1} = ${D(S.Mxp, 1000, 1)} kNm/m`);
      R('Moment capacity of plate', `M_{c} = 1.2 × Z_{x,plate} × p_{ysb} = ${D(S.Mc, 1000, 1)} kNm/m`);
      C(ok(u.heel) ? 'PASS - Design strength exceeds stress at heel' : 'FAIL - Stress at heel exceeds design strength');   // FAIL: live reference run V2
      const s1 = S.Mx / S.Zxx, cfp = (4 * (o.pysb * 1e6) ** 2 - 3 * s1 ** 2) ** 0.5, ats = Math.max((cfp ** 2 - s1 ** 2) / (2 * cfp * o.pysb * 1e6), 0);
      R('Longitudinal stress due to overall bending', `σ_{1} = M_{x} / Z_{xx,all} = ${D(s1, 1e6, 1)}${N2}`);
      R('Constant relating to Von Mises curve', `c_{fp} = (4 × p_{ysb}^{2} - 3 × σ_{1}^{2})^{0.5} = ${D(cfp, 1e6, 1)}${N2}`);
      R('Transverse bending stress ratio limit', `α_{ts} = (c_{fp}^{2} - σ_{1}^{2}) / (2 × c_{fp} × p_{ysb}) = ${D(ats, 1, 3)}`);
      R('Transverse bending stress ratio', `α_{ls} = M_{x,plate} / M_{c} = ${D(S.als, 1, 3)}`);
      C(ok(u.transverse) ? 'PASS - Transverse bending stress ratio less than allowable limit' : 'FAIL - Transverse bending stress ratio exceeds allowable limit');
      H('Deflection at toe');
      R('Unfactored force on support angle', `${P1Sf} = ${D(S.P1SLS, 1000, 1)} kN/m`);
      R('Distance from weld to load position', `a_{m} = e_{1} = ${D(S.am, 1e-3, 0)} mm`);
      R('Length of load resultant to edge of plate', `b_{m} = l_{h} - e_{1} = ${D(S.bm, 1e-3, 0)} mm`);
      R('Dist from weld to load position as ratio of length', `a_{l} = a_{m} / (a_{m} + b_{m}) = ${D(S.al, 1, 3)}`);
      R('Effective second moment of inertia', `I_{eff_def} = t_{sb}^{3} / 12 = ${D(S.Ieff, 1e-12, 0)} mm^{4}/m`);
      R('Deflection at toe', `δ = (a_{l}^{2} × (3 - a_{l}) / 6) × (P_{1SLS} × (a_{m} + b_{m})^{3}) / (E_{S5950} × I_{eff_def}) = ${D(S.delta, 1e-3, 2)} mm`);
      R('Deflection limit', `δ_{lim} = ${D(S.dlim, 1e-3, 2)} mm`);
      C(ok(u.toe) ? 'PASS - Deflection is within specified criteria' : 'FAIL - Deflection exceeds specified criteria');
      if (p.con_stage) {
        H('Construction stage biaxial stress effects in the plate');
        R('Eccentricity of support beam masonry', `e_{1c} = c - e_{mb} + b_{msb} / 2 - s = ${D(S.e1c, 1e-3, 0)} mm`);
        R('Force of masonry on support plate', `P_{1c} = 1.2m × b_{msb} × ρ_{m,sb} × γ_{fG} = ${D(S.P1c, 1000, 1)} kN/m`);
        R('Bending at heel', `M_{x,platec} = P_{1c} × e_{1c} = ${D(S.Mxpc, 1000, 1)} kNm/m`);
        C(ok(u.heel_c) ? 'PASS - Design strength exceeds stress at heel' : 'FAIL - Stress at heel exceeds design strength');
        R('Transverse bending stress ratio', `α_{lsc} = M_{x,platec} / M_{c} = ${D(S.alsc, 1, 3)}`);
        C(ok(u.transverse_c) ? 'PASS - Transverse bending stress ratio less than allowable limit' : 'FAIL - Transverse bending stress ratio exceeds allowable limit');
        H('Construction stage deflection at toe');
        R('Unfactored force on support angle', `P_{1cSLS} = 1.2m × b_{msb} × ρ_{m,sb} = ${D(S.P1cS, 1000, 1)} kN/m`);
        R('Distance from weld to load position', `a_{mc} = e_{1c} = ${D(S.amc, 1e-3, 0)} mm`);
        R('Length of load resultant to edge of plate', `b_{mc} = l_{h} - e_{1c} = ${D(S.bmc, 1e-3, 0)} mm`);
        R('Dist from weld to load position as ratio of length', `a_{lc} = a_{mc} / (a_{mc} + b_{mc}) = ${D(S.alc, 1, 3)}`);
        R('Deflection at toe', `δ_{c} = (a_{lc}^{2} × (3 - a_{lc}) / 6) × (P_{1cSLS} × (a_{mc} + b_{mc})^{3}) / (E_{S5950} × I_{eff_def}) = ${D(S.deltac, 1e-3, 2)} mm`);
        C(ok(u.toe_c) ? 'PASS - Deflection is within specified criteria' : 'FAIL - Deflection exceeds specified criteria');
      }
      H('Weld details - assume a full length weld and that the plate acts as a propped cantilever with the prop at the weld position and the fixed end at the centre of the torsion beam');
      R('Leg length of weld', `s_{weld} = ${MM(p.s_weld)} mm`);
      R('Throat size of weld', `a_{weld} = 1/√(2) × s_{weld} = ${D(S.aw, 1e-3, 1)} mm`);
      R('Shear force at weld position', `R_{A} = P_{1} × maX((1 + (3 × e_{1}) / (2 × B_{mb} / 2)), 1.4) = ${D(S.RA, 1000, 1)} kN/m`);
      R('Maximum possible force in plate', `R_{p} = (l_{h} + B_{mb}) × t_{sb} × p_{ysb} = ${D(S.Rp, 1000, 1)} kN`);
      R('Longitudinal shear between beam and plate', `R_{l} = 2 × R_{p} / L = ${D(S.Rl, 1000, 1)} kN/m`);
      R('Horizontal shear between beam and plate', `R_{h} = P_{1} × e_{1} / (s_{weld} / 2 + t_{sb} / 2) = ${D(S.Rh, 1000, 1)} kN/m`);
      R('Resultant weld force', `R_{weld} = (R_{A}^{2} + R_{l}^{2} + R_{h}^{2})^{0.5} = ${D(S.Rweld, 1e6, 3)} kN/mm`);
      R('Strength of weld (Table 37)', `p_{weld} = ${D(S.pweld, 1e6, 1)}${N2}`);
      R('Capacity of full length weld', `p_{c,weld} = a_{weld} × p_{weld} = ${D(S.pcw, 1e6, 3)} kN/mm`);
      C(ok(u.weld) ? 'PASS - Capacity of weld exceeds resultant force on weld' : 'FAIL - Resultant force on weld exceeds capacity of weld');   // FAIL: live reference run V7
    } else if (m === 'Cantilever') {
      // ---- hot-rolled angle (SCI P157)
      H('Forces simple');
      const where = bolt ? 'support angle' : 'support beam';
      R(`Eccentricity of ${where} ${what}`, `e_{1} = ${D(S.e1, 1e-3, 0)} mm`);
      R(`Force of ${what} on ${where}`, `${P1f} = ${D(S.P1, 1000, 1)} kN/m`);
      R('Bending at heel', `M_{1} = P_{1} × e_{1} = ${D(S.M1, 1000, 1)} kNm/m`);
      R('Section modulus of horizontal plate', `Z_{xx,plate} = 1m × t_{sb}^{2} / (6 × 1m) = ${D(S.Zp, 1e-6, 2)} cm^{3}/m`);
      R('Moment capacity of plate', `M_{c} = 1.2 × p_{ysb} × Z_{xx,plate} = ${D(S.Mc, 1000, 1)} kNm/m`);
      if (bolt) {
        R('Effective width factor', `λ_{m} = B_{bolt} / (l_{h} + l_{b}) = ${D(S.lamm, 1, 2)}`);
        R('Effect. breadth of angle between fixings for stress', `B_{eff_str} = min(2 × B_{bolt} / (0.4 + λ_{m}), B_{bolt}) = ${D(S.Bstr, 1e-3, 0)} mm`);
        R('Stress at heel', `σ = B_{bolt} × M_{1} / (B_{eff_str} × Z_{xx,plate}) = ${D(S.sigma, 1e6, 1)}${N2}`);
      } else R('Stress at heel', `σ = M_{1} / Z_{xx,plate} = ${D(S.sigma, 1e6, 1)}${N2}`);
      C(ok(u.heel) ? 'PASS - Design strength exceeds stress at heel' : 'FAIL - Stress at heel exceeds design strength');
      H('Deflection at toe');
      R('Unfactored force on support angle', `${P1Sf} = ${D(S.P1SLS, 1000, 1)} kN/m`);
      R('Circumferential distance bearing point to load pos', `a_{m} = e_{1}${hb < 0 ? ' + h_{bf}' : ''} = ${D(S.am, 1e-3, 0)} mm`);
      R('Length of load resultant to toe edge of angle', `b_{m} = l_{h} - e_{1} = ${D(S.bm, 1e-3, 0)} mm`);
      R('Ratio p', `p_{def} = a_{m} / l = ${D(S.pdef, 1, 2)}`);
      R('Ratio q', `q_{def} = b_{m} / l = ${D(S.qdef, 1, 2)}`);
      if (bolt) {
        items.push({ t: 'row', l: 'Effective breadth of angle between fixings for defl.', e: `B_{eff_def} = min(2 × B_{bolt} / (0.8 + λ_{m}), B_{bolt}) = ${D(S.Bdef, 1e-3, 0)} mm` });
        R('Effective second moment of inertia', `I_{eff_def} = B_{eff_def} × t_{sb}^{3} / 12 = ${D(S.Ieff, 1e-12, 0)} mm^{4}`);
        R('Deflection at toe', `δ = (B_{bolt} × P_{1SLS} × l^{3} × p_{def}) / (12 × E_{S5950} × I_{eff_def}) × (4 × p_{def}^{2} + 6 × p_{def} × q_{def} + 3 × p_{def} + 3 × q_{def}) = ${D(S.delta, 1e-3, 2)} mm`);
      } else {
        R('Effective second moment of inertia', `I_{eff_def} = 1m × t_{sb}^{3} / (12 × 1m) = ${D(S.Ieff, 1e-12, 0)} mm^{4}/m`);
        R('Deflection at toe', `δ = (P_{1SLS} × l^{3} × p_{def}) / (12 × E_{S5950} × I_{eff_def}) × (4 × p_{def}^{2} + 6 × p_{def} × q_{def} + 3 × p_{def} + 3 × q_{def}) = ${D(S.delta, 1e-3, 2)} mm`);
      }
      R('Deflection limit', `δ_{lim} = ${D(S.dlim, 1e-3, 2)} mm`);
      C(ok(u.toe) ? 'PASS - Deflection is within specified criteria' : 'FAIL - Deflection exceeds specified criteria');
      if (!bolt) {
        H('Weld details');
        R('Leg length of weld', `s_{weld} = ${MM(p.s_weld)} mm`);
        R('Throat size of weld', `a_{weld} = 1/√(2) × s_{weld} = ${D(S.aw, 1e-3, 1)} mm`);
        R('Strength of weld (Table 37)', `p_{weld} = ${D(o.pweld * 1e6, 1e6, 0)}${N2}`);
        R('Resultant force on weld', `f_{weld} = √(P_{1}^{2} + (M_{1} / l)^{2}) = ${D(S.fweld, 1000, 1)} kN/m`);
        R('Provided length of weld', `l_{weld} = ${B(p.l_weld * 1000, 0)} mm/m`);
        R('Resultant stress on weld', `P_{weld} = f_{weld} / (l_{weld} × a_{weld}) = ${D(S.Pweld, 1e6, 0)}${N2}`);
        C(ok(u.weld) ? 'PASS - Strength of weld exceeds resultant stress on weld' : 'FAIL - Resultant stress on weld exceeds strength of weld');   // FAIL: live reference run V3
      }
    } else {
      // ---- cold-formed stainless angle (SCI P157)
      H('Stress in heel');
      R('Force of masonry on support beam', `${P1f} = ${D(S.P1, 1000, 1)} kN/m`);
      R('Length from bolt CL to start of angle radius', `l = l_{b} - (t_{sb} + r_{1sb}) = ${D(S.l, 1e-3, 0)} mm`);
      R('Bolt edge distance', `g = l_{v} - l_{b} = ${D(S.g, 1e-3, 0)} mm`);
      R('Position of resultant load', `e = max(10mm, min(b_{msb} / 2, min(((16.5 × t_{sb} / (l + g)) - 0.73), 1) × d_{m})) = ${D(S.e, 1e-3, 0)} mm`);
      R('Circumferential distance bearing point to load pos', `a_{m} = c${tor ? ' - e_{mb}' : ''} - (s + t_{sb} + r_{1sb}) + π × (t_{sb} / 2 + r_{1sb}) / 2 + e = ${D(S.am, 1e-3, 0)} mm`);
      R('Length of load resultant to toe edge of angle', `b_{m} = d_{m} - e = ${D(S.bm, 1e-3, 0)} mm`);
      R('Effective width factor', `λ_{m} = B_{bolt} / (l + a_{m} + b_{m}) = ${D(S.lamm, 1, 2)}`);
      R('Effect. breadth of angle between fixings for stress', `B_{eff_str} = min(2 × B_{bolt} / (0.4 + λ_{m}), B_{bolt}) = ${D(S.Bstr, 1e-3, 0)} mm`);
      R('Stress at heel of angle', `σ = 6 × P_{1} × B_{bolt} × (c${tor ? ' - e_{mb}' : ''} + e - s - 0.5 × t_{sb}) / (B_{eff_str} × t_{sb}^{2}) = ${D(S.sigma, 1e6, 1)}${N2}`);
      C(ok(u.heel) ? 'PASS - Proof stress exceeds stress at heel' : 'FAIL - Stress at heel exceeds proof stress');   // FAIL: live reference run V4
      H('Deflection at toe');
      R('Unfactored force on support angle', `${P1Sf} = ${D(S.P1SLS, 1000, 1)} kN/m`);
      items.push({ t: 'row', l: 'Effective breadth of angle between fixings for defl.', e: `B_{eff_def} = min(2 × B_{bolt} / (0.8 + λ_{m}), B_{bolt}) = ${D(S.Bdef, 1e-3, 0)} mm` });
      R('Effective second moment of inertia', `I_{eff_def} = B_{eff_def} × t_{sb}^{3} / 12 = ${D(S.Ieff, 1e-12, 0)} mm^{4}`);
      R('Ratio p', `p_{def} = a_{m} / l = ${D(S.pdef, 1, 2)}`);
      R('Ratio q', `q_{def} = b_{m} / l = ${D(S.qdef, 1, 2)}`);
      R('Deflection at toe', `δ = (B_{bolt} × P_{1SLS} × l^{3} × p_{def}) / (12 × E_{S5950} × I_{eff_def}) × (4 × p_{def}^{2} + 6 × p_{def} × q_{def} + 3 × p_{def} + 3 × q_{def}) = ${D(S.delta, 1e-3, 2)} mm`);
      R('Deflection limit', `δ_{lim} = min((1 + d_{m} / c) × 1mm, 2mm) = ${D(S.dlim, 1e-3, 2)} mm`);
      C(ok(u.toe) ? 'PASS - Deflection is within specified criteria' : 'FAIL - Deflection exceeds specified criteria');
    }
    if (bolt) {
      H('Bolt details');
      R('Type of bolt', `Grade ${p.bolt_grade} Black Bolts - ${p.bolt_size}`);
      R('Effective area of bolt', `A_{sbolt} = ${B(o.As_bolt, 1)} mm^{2}`);
      R('Shear capacity of bolt', `P_{sbolt} = ${B(o.Ps_bolt / 1e3, 1)} kN`);
      R('Shear force on bolt', `F_{s,bolt} = P_{1} × B_{bolt} = ${D(S.Fs, 1000, 3)}kN`);
      R('Tension capacity of bolt', `P_{tbolt} = ${B(o.Pt_bolt / 1e3, 1)} kN`);
      R('Tension force on bolt', `F_{t,bolt} = ${m === 'SCI' ? 'B_{bolt} × P_{1} × a_{m} / l' : 'B_{bolt} × M_{1} / l'} = ${D(S.Ft, 1000, 3)}kN`);
      C(ok(u.bolt) ? 'PASS - Capacity of bolt exceeds applied force' : 'FAIL - Applied force exceeds capacity of bolt');   // FAIL: live reference run V4
      R('Combined shear and tension', `F_{s,bolt} / P_{sbolt} + F_{t,bolt} / P_{tbolt} = ${D(S.Fcomb, 1, 3)}`);
      C(ok(u.bolt_comb) ? 'PASS - Combined shear and tension check is less than 1.4' : 'FAIL - Combined shear and tension check is greater than 1.4');   // live reference run V2
    }
    if (!tor) return items;

    // ---- torsional loading of the beam
    const isC = mb.type === 'C', hr = mb.type === 'HR';
    const A3 = m === 'Cantilever' ? 'A_{sb}' : 'A_{sbu}';
    H('Torsional loading ULS');
    R(dirS ? 'Loading of support beam' : 'Loading of support beam masonry', `w_{1ULS} = ${slabS ? 'P_{Gssb} × γ_{fG} + P_{Qssb} × γ_{fQ}' : dirS ? 'G_{sb} × γ_{fG} + Q_{sb} × γ_{fQ}' : '(h_{msb} × b_{msb} × ρ_{m,sb} + P_{Gaddsb}) × γ_{fG} + P_{Qaddsb} × γ_{fQ}'} = ${D(S.w1, 1000, 2)} kN/m`);
    R(dirM ? 'Loading of main beam' : 'Loading of main beam masonry', `w_{2ULS} = ${slabM ? 'P_{Gsmb} × γ_{fG} + P_{Qsmb} × γ_{fQ}' : dirM ? 'G_{mb} × γ_{fG} + Q_{mb} × γ_{fQ}' : '(h_{mmb} × b_{mmb} × ρ_{m,mb} + P_{Gaddmb}) × γ_{fG} + P_{Qaddmb} × γ_{fQ}'} = ${D(S.w2, 1000, 2)} kN/m`);
    R('Self weight of support beam', `w_{3ULS} = ${A3} × ρ_{sb} × γ_{fG} = ${D(S.w3, 1000, 2)} kN/m`);
    H('Torsional loading SLS');
    R(dirS ? 'Loading of support beam' : 'Loading of support beam masonry', `w_{1SLS} = ${slabS ? 'P_{Gssb} + P_{Qssb}' : dirS ? 'G_{sb} + Q_{sb}' : 'h_{msb} × b_{msb} × ρ_{m,sb} + P_{Gaddsb} + P_{Qaddsb}'} = ${D(S.w1S, 1000, 2)} kN/m`);
    R(dirM ? 'Loading of main beam' : 'Loading of main beam masonry', `w_{2SLS} = ${slabM ? 'P_{Gsmb} + P_{Qsmb}' : dirM ? 'G_{mb} + Q_{mb}' : 'h_{mmb} × b_{mmb} × ρ_{m,mb} + P_{Gaddmb} + P_{Qaddmb}'} = ${D(S.w2S, 1000, 2)} kN/m`);
    R('Self weight of support beam', `w_{3SLS} = ${A3} × ρ_{sb} = ${D(S.w3S, 1000, 2)} kN/m`);
    H('Eccentricities');
    R('Distance to shear centre of main beam', `e_{0mb} = ${B(o.e0mb, 0)} mm`);
    const face = isC ? 't_{mb} / 2' : 'B_{mb} / 2';
    if (slabS) R('Eccentricity of support beam slab', `e_{1mb} = ${isC ? 'e_{1} + s + t_{mb} / 2 - e_{0mb}' : 'e_{1} - e_{0mb} + s + B_{mb} / 2'} = ${D(S.e1mb, 1e-3, 0)} mm`);   // channel: live reference run V8
    else R(`Eccentricity of support beam ${dirS ? 'load' : 'masonry'}`, `e_{1mb} = ${isC ? 'c + (b_{msb} + t_{mb}) / 2 - e_{mb} - e_{0mb}' : '(B_{mb} + b_{msb}) / 2 + c - e_{mb}'} = ${D(S.e1mb, 1e-3, 0)} mm`);
    if (slabM) R('Eccentricity of main beam slab', `e_{2mb} = ${isC ? '- e_{2} - e_{0mb} + t_{mb} / 2' : 'B_{mb} / 2 - e_{2}'} = ${D(S.e2mb, 1e-3, 0)} mm`);   // channel: live reference run V8
    else R(`Eccentricity of main beam ${dirM ? 'load' : 'masonry'}`, `e_{2mb} = ${isC ? '-(b_{mmb} - t_{mb}) / 2 - e_{mb} - e_{0mb}' : '(B_{mb} - b_{mmb}) / 2 - e_{mb}'} = ${D(S.e2mb, 1e-3, 0)} mm`);
    R('Eccentricity of support beam', `e_{3mb} = ${isC ? 'c_{yysb} - e_{0mb} + t_{mb} / 2' : 'B_{mb} / 2 + c_{yysb}'} = ${D(S.e3mb, 1e-3, 0)} mm`);
    H('Torsional effects');
    R('Applied torque (ULS)', `T_{qULS} = abs(w_{1ULS} × e_{1mb} + w_{2ULS} × e_{2mb} + w_{3ULS} × e_{3mb}) = ${D(S.TqULS, 1000, 2)} kNm/m`);
    R('Total torque (ULS)', `T_{q} = T_{qULS} × L = ${D(S.Tq, 1000, 2)} kNm`);
    R('Applied torque (SLS)', `T_{qSLS} = abs(w_{1SLS} × e_{1mb} + w_{2SLS} × e_{2mb} + w_{3SLS} × e_{3mb}) = ${D(S.TqSLS, 1000, 2)} kNm/m`);
    R('Total torque (SLS)', `T_{qu} = T_{qSLS} × L = ${D(S.Tqu, 1000, 2)} kNm`);

    // ---- steel beam torsion design
    items.push({ t: 'title', s: 'STEEL BEAM TORSION DESIGN' });
    items.push({ t: 'sub', s: 'In accordance with BS5950-1:2000 incorporating Corrigendum No.1' });
    items.push({ t: 'ver', s: opts.versionLine2 || '' });
    H('Section details');
    R('Section type', mb.name);
    R('Steel grade', p.grade_mb);
    R('Design stength', `p_{yw} = p_{y} = ${B(o.py, 0)}${N2}`);
    R('Constant', `ε = √(275 N/mm^{2} / p_{y}) = ${D(S.eps, 1, 3)} `);
    H('Geometry - Beam unrestrained against lateral-torsional buckling between supports. ');
    R('Effective span', `L = ${MM(p.L)} mm`);
    R('Length of segment for LT buckling', `L_{LT} = ${MM(p.LLT)} mm`);
    for (const s of RESTRAINT[restraintOf(p)]) T(s);
    R('Effective length for LT buckling', `L_{E_LT} = L_{LT} × ${kText(p.K_LT)} = ${D(S.LE, 1e-3, 0)} mm`);
    H('Loading - Torsional loading comprises only full-length uniformly distributed load(s)');
    H('Internal forces & moments on member under factored loading for uls design');
    R('Applied shear force', `F_{vy} = ${D(S.Fv, 1000, 1)} kN`);
    R('Maximum bending moment', `M_{LT} = M_{x} = ${D(S.Mx, 1000, 2)} kNm`);
    R('Applied torque', `T_{q} = ${D(S.Tq, 1000, 2)} kNm`);
    R('Minor axis bending moment', 'M_{y} = 0 kNm');            // plain in the reference
    R('Compression force', 'F_{c} = 0 kN');
    H('Equivalent uniform moment factors');
    R('EUM factor (Cl. 4.3.6.6 and T18)', `m_{LT} = ${B(p.mLT, 3)}`);
    const DEG = 180 / PI_T;
    if (hr) {
      H('Torsional deflection analysis');
      T('Beam is torsion fixed at each end. (as defined in SCI-P-057 section 2.1.6)');
      R('Maximum torque (at supports)', `T_{o} = T_{q} / 2 = ${D(S.Tq / 2, 1000, 2)} kNm`);
      R('Average torque between support & centreline', `T_{av} = T_{o} / 2 = ${D(S.Tq / 2 / 2, 1000, 2)} kNm`);
      R('Max. angle of twist (at midspan)', `φ = T_{av} / (G × J) × L / 2 × 1 rads = ${B(o.phi_deg / DEG, 3)} rads`);
    } else {
      const sg = o.twistSigned || {};
      H('Torsional deflection parameters');
      if (p.warping === 'Free') {
        T('Beam is torsion fixed and warping free at each end. (as defined in SCI-P-057 section 2.1.6) - Appendix B case 4');
        R('Dist along the beam for first derivative of twist', 'z_{1} = 0 mm');
        R2('Dist along the beam for second derivative of twist', `z_{2} = L / 2 = ${B(p.L / 2, 0)} mm`);
        R2('First derivative of angle of twist', `φ'_{1} = T_{q} / (G × J) × a / L × [L^{2} / (2 × a) × (1 / L - 2 × z_{1} / L^{2}) + sinh(z_{1} / a) - tanh(L / (2 × a)) × cosh(z_{1} / a)] × 1 rads = ${sci(sg.p1 * 1e3)} rads/m`);
        R2('Third derivative of angle of twist', `φ'''_{1} = T_{q} / (G × J × a^{2}) × a/L × [sinh(z_{1} / a) - tanh(L / (2 × a)) × cosh(z_{1} / a)] × 1 rads = ${sci(sg.p3 * 1e9)} rads/m^{3}`);
        R2('Angle of twist', `φ_{2} = T_{q} × a / (G × J) × a / L × [L^{2} / (2 × a^{2}) × (z_{2} / L - z_{2}^{2} / L^{2}) + cosh(z_{2} / a) - tanh(L / (2 × a)) × sinh(z_{2} / a) - 1] × 1 rads = ${B(sg.phi, 3)} rads`);
        R2('Second derivative of angle of twist', `φ''_{2} = T_{q} / (G × J × a) × a / L × [cosh(z_{2} / a) - tanh(L / (2 × a)) × sinh(z_{2} / a) - 1] × 1 rads = ${sci(sg.p2 * 1e6)} rads/m^{2}`);
        H('Design parameters');
        R('Total angle of twist', `φ = abs(φ_{2}) = ${B(Math.abs(sg.phi), 3)} rads`);
        R('First derivative of φ', `φ' = abs(φ'_{1}) = ${sci(Math.abs(sg.p1) * 1e3)} rads/m`);
        R('Second derivative of φ', `φ'' = abs(φ''_{2}) = ${sci(Math.abs(sg.p2) * 1e6)} rads/m^{2}`);
        R('Third derivative of φ', `φ''' = abs(φ'''_{1}) = ${sci(Math.abs(sg.p3) * 1e9)} rads/m^{3}`);
      } else {
        T('Beam is torsion fixed and warping fixed at each end. (as defined in SCI-P-057 section 2.1.6) - Appendix B case 6');
        R('Dist along the beam for first derivative of twist', 'z_{1} = 0 mm');
        R('Dist along the beam for second derivative of twist', `z_{2} = 0.2 × L = ${B(0.2 * p.L, 0)} mm`);
        R('Dist along the beam for third derivative of twist', `z_{3} = L / 2 = ${B(p.L / 2, 0)} mm`);
        R('Second derivative of angle of twist', `φ''_{1} = T_{q} / (2 × G × J × a) × [ (1 + cosh(L / a)) / sinh(L / a) × cosh(z_{1} / a) - 2 × a / L - sinh(z_{1} / a)] × 1 rads = ${sci(sg.p2 * 1e6)} rads/m^{2}`);
        R('Third derivative of angle of twist', `φ'''_{1} = T_{q} / (2 × G × J × a^{2}) × [(1 + cosh(L / a)) / sinh(L / a) × sinh(z_{1} / a) - cosh(z_{1} / a)] × 1 rads = ${sci(sg.p3 * 1e9)} rads/m^{3}`);
        R('First derivative of angle of twist', `φ'_{2} = T_{q} / (2 × G × J) × [(1 + cosh(L / a)) / sinh(L / a) × sinh(z_{2} / a) + (1 - 2 × z_{2} / L) - cosh(z_{2} / a)] × 1 rads = ${sci(sg.p1 * 1e3)} rads/m`);
        R('Angle of twist', `φ_{3} = T_{q} × a / (2 × G × J) × [(1 + cosh(L / a)) / sinh(L / a) × (cosh(z_{3} / a) - 1) + z_{3} / a × (1 - z_{3} / L) - sinh(z_{3} / a)] × 1 rads = ${B(sg.phi, 3)} rads`);
        H('Design parameters');
        R('Total angle of twist', `φ = abs(φ_{3}) = ${B(Math.abs(sg.phi), 3)} rads`);
        R('First derivative of φ', `φ' = abs(φ'_{2}) = ${sci(Math.abs(sg.p1) * 1e3)} rads/m`);
        R('Second derivative of φ', `φ'' = abs(φ''_{1}) = ${sci(Math.abs(sg.p2) * 1e6)} rads/m^{2}`);
        R('Third derivative of φ', `φ''' = abs(φ'''_{1}) = ${sci(Math.abs(sg.p3) * 1e9)} rads/m^{3}`);
      }
    }
    H('Section classification');
    if (hr) {
      R('', `b_{x} / t = ${B(o.bt, 1)}`);
      R('', `d_{x} / t = ${B(o.dt, 1)}`);
      R('', `b_{y} / t = ${B(o.dt, 1)}`);
      R('', `d_{y} / t = ${B(o.bt, 1)}`);
      R('', `r_{1sx} = min(1.0, max(-1.0, F_{c} / ( 2 × d_{x}× t × p_{yw}))) = ${B(0, 3)}`);
      R('', `r_{1sy} = min(1.0, max( -1.0, F_{c} / ( 2 × d_{y}× t × p_{yw}))) = ${B(0, 3)}`);
      R('', `r_{2s} = F_{c} / (A_{g} × p_{yw}) = ${B(0, 3)}`);
    } else {
      R('', `b / T = ${B(o.bt, 1)}`);
      R('', `d / t = ${B(o.dt, 1)}`);
      R('', `r_{1s} = min( 1.0, max(-1.0, F_{c} / (d × t × p_{yw} ))) = ${B(0, 3)}`);
      R('', `r_{2s} = F_{c} / (A_{g} × p_{yw}) = ${B(0, 3)}`);
    }
    C(`Section classification is ${String(o.Class).toLowerCase()}`);
    if (o.slender) {
      C('Section is slender - beyond scope');   // UNVERIFIED placement
      return items;
    }
    const x = o.util;
    H('Shear capacity (parallel to y-axis) ');
    R('Design shear force', `F_{vy} = ${D(S.Fv, 1000, 1)} kN`);
    R('Design shear resistance (Cl. 4.2.3)', `P_{vy} = 0.6 ×  p_{y} ×  A_{vy} = ${D(S.Pv, 1000, 1)} kN`);
    C(ok(x.shear) ? 'Pass - Shear' : 'FAIL - Shear');   // FAIL: live reference run V5
    const plastic = o.Class === 'Plastic' || o.Class === 'Compact';
    H('Moment capacity (x-axis) ');
    R('Design bending moment', `M_{x} = ${D(S.Mx, 1000, 1)} kNm`);
    if (o.highShear) {   // live reference runs V1 (UKC) and V3 (RHS), plastic; semi-compact not sampled
      R('Reduction factor', `ρ_{y} = (2 × (F_{vy} / P_{vy}) -1)^{2} = ${D(S.rhoy, 1, 2)}`);
      R('Moment capacity', `M_{cxu} = p_{y} × (S_{x} - ρ_{y} × S_{vx}) = ${D(S.Mcxu, 1000, 1)} kNm`);
      R('Moment capacity high shear (Cl. 4.2.5.1)', `M_{cx} = min(p_{y} × (S_{x} - ρ_{y} × S_{vx}), 1.2 × p_{y} × Z_{x}) = ${D(S.Mcx, 1000, 1)} kNm`);
    } else {
      R('Moment capacity', `M_{cxu} = p_{y} × ${plastic ? 'S' : 'Z'}_{x} = ${D(S.Mcxu, 1000, 1)} kNm`);
      R('Moment capacity low shear (Cl. 4.2.5.1)', `M_{cx} = min(p_{y} × ${plastic ? 'S' : 'Z'}_{x}, 1.2 × p_{y} × Z_{x}) = ${D(S.Mcx, 1000, 1)} kNm`);
    }
    C(ok(x.moment) ? 'Pass - Moment capacity exceeds design bending moment' : 'FAIL - Design bending moment exceeds moment capacity');
    H('Lateral torsional buckling ');
    if (hr) {
      if (mb.D <= mb.B) {
        T('LT buckling check not required for this section (cl. 4.6.3.1)');   // style UNVERIFIED (no PDF sample)
      } else {
        R('Effective length for lateral torsional buckling', `L_{E_LT} = ${D(S.LE, 1e-3, 0)} mm`);
        R('Slenderness ratio - cl 4.3.6.5', `λ = L_{E_LT} / r_{y} = ${D(S.lam, 1, 0)}`);
        R('', `D / B = ${B(mb.D / mb.B, 1)}`);
        C(o.phib === undefined ? 'LTB check not required' : 'LTB check is required');
      }
      if (o.phib !== undefined) {   // beyond the Table 15 limit: Annex B.2.6.1, as the reference prints it (live run V6)
        R('', `γ_{b} = (1 - I_{y} / I_{x}) × (1 - J / (2.6 × I_{x})) = ${D(S.gb, 1, 2)}`);
        R('', `φ_{b} =  √(S_{x}^{2} × γ_{b} / (A × J)) = ${D(S.phib, 1, 2)}`);
        R('Ratio - cl 4.3.6.9', plastic ? `β_{w} =  1.0 = ${B(1, 3)}` : `β_{w} =  Z_{x} / S_{x} = ${D(S.bw, 1, 3)}`);   // semi-compact form UNVERIFIED
        R('Equvalent slenderness – cl 4.3.6.7', `λ_{LT} =  2.25 × √(φ_{b} × λ × β_{w}) = ${D(S.lamLT, 1, 0)}`);
        R('Limiting slendernes – Annex B2.2', `λ_{L0} =  0.4 × √(π^{2} × E_{S5950} / p_{y}) = ${D(S.lamL0, 1, 0)}`);
        R('Euler stress', `p_{E} =  π^{2} × E_{S5950} / λ_{LT}^{2} = ${D(S.pE, 1e6, 0)}${N2}`);
        R('Perry factor', `η_{LT} =  max(7.0 × (λ_{LT} - λ_{L0}) / 1000, 0) = ${D(S.eta, 1, 3)}`);
        R('', `φ_{LT} =  (p_{y} + (η_{LT} + 1 ) × p_{E} ) / 2 = ${D(S.phiLT, 1e6, 0)}${N2}`);
        R('Bending strength', `p_{b} =  p_{E} × p_{y} / (φ_{LT} + √(φ_{LT}^{2} - p_{E}× p_{y})) = ${D(S.pb, 1e6, 0)}${N2}`);
        R('Buckling resistance moment', `M_{b} =  p_{b} × ${plastic ? 'S' : 'Z'}_{x} = ${D(S.Mb, 1000, 1)} kNm`);
        R('Max moment governing buckling resistance', `M_{LT} = ${D(S.Mx, 1000, 1)} kNm`);
        R('Equiv uniform moment factor for LTB', `m_{LT} = ${B(p.mLT, 2)}`);
        R('', `M_{b} / m_{LT} = ${D(S.Mb / p.mLT, 1000, 1)} kNm`);
        C(ok(x.LTB) ? 'Pass - lat. tors. buckling' : 'FAIL - lat. tors. buckling');
      } else {
        R('Buckling resistance moment', `M_{b} = M_{cx} = ${D(S.Mb, 1000, 1)} kNm`);
        if (mb.D <= mb.B) C('LT buckling check not required for this section');
      }
    } else {
      R('Effective length for lateral torsional buckling', `L_{E_LT} = ${D(S.LE, 1e-3, 0)} mm`);
      R('Slenderness ratio', `λ = L_{E_LT} / r_{y} =${D(S.lam, 1, 0)}`);
      R('Buckling parameter', `u = ${B(mb.u, 3)}`);
      R('Flange ratio', 'η = 0.5');
      R('Torsional index', `x = ${B(mb.x, 1)}`);
      R('Slenderness factor', `v = 1 / (1 + 0.05 × (λ / x)^{2})^{0.25} = ${D(S.v, 1, 2)}`);
      R('Ratio - cl 4.3.6.9', plastic ? `β_{w} = 1.0 = ${B(1, 3)}` : `β_{w} = Z_{x} / S_{x} = ${D(S.bw, 1, 3)}`);
      R('Equvalent slenderness – cl 4.3.6.7', `λ_{LT} = u × v × λ × √(β_{w}) = ${D(S.lamLT, 1, 0)}`);
      R('Limiting slendernes – Annex B2.2', `λ_{L0} = 0.4 × √(π^{2}× E_{S5950} / p_{y}) = ${D(S.lamL0, 1, 0)}`);
      R('Euler stress', `p_{E} = π^{2} × E_{S5950} / λ_{LT}^{2} = ${D(S.pE, 1e6, 0)}${N2}`);
      R('Perry factor', `η_{LT} = max(7.0 × ( λ_{LT} - λ_{L0}) / 1000, 0) = ${D(S.eta, 1, 3)}`);
      R('', `φ_{LT} = (p_{y} + (η_{LT} + 1) × p_{E}) / 2 = ${D(S.phiLT, 1, 3)}`);
      R('Bending strength', `p_{b} = p_{E} × p_{y} / (φ_{LT} + √(φ_{LT}^{2} - p_{E }× p_{y})) = ${D(S.pb, 1e6, 0)}${N2}`);   // 'E ' subscripted, as the reference
      R('Buckling resistance moment', `M_{b} = p_{b} × ${plastic ? 'S' : 'Z'}_{x} = ${D(S.Mb, 1000, 1)} kNm`);
      R('Max moment governing buckling resistance', `M_{LT} = ${D(S.Mx, 1000, 1)} kNm`);
      R('Equiv uniform moment factor for LTB', `m_{LT} = ${B(p.mLT, 2)}`);
      R('', `M_{b} / m_{LT} = ${D(S.Mb / p.mLT, 1000, 1)} kNm`);
      C(ok(x.LTB) ? 'Pass - lat. tors. buckling' : 'FAIL - lat. tors. buckling');
    }
    const phi = o.phi_deg / DEG;
    if (hr) {
      H('Buckling under combined bending & torsion - SCI-P-057 section 2.3');
      T('For simplicity, a conservative check is applied using the maximum stresses due to each of the separate load effects, even though these do not necessarily all occur at the same section along the member.');
      R('Maximum angle of twist', `φ = ${B(phi, 3)} rads`);
      R('Induced minor axis moment', `M_{yt} = M_{x} × φ / 1 rad= ${D(S.Myt, 1000, 2)} kNm`);
      R('Normal stress at corner due to M_{yt}', `σ_{byt} = M_{yt} / Z_{y} = ${D(S.sbyt, 1e6, 0)}${N2}`);
      R('Interaction index', `i_{b} = M_{x} × m_{LT} / M_{b} + σ_{byt} / p_{y} × ( 1 + 0.5 × M_{x} × m_{LT} / M_{b}) = ${B(o.ib, 2)}`);
    } else {
      H('Buckling under combined bending & torsion -SCI-P-057 section 2.3');
      T('For simplicity, a conservative check is applied using the maximum stresses due to each of the separate load effects, even though these do not necessarily all occur at the same section along the member.');
      R('Span factor', `L / a = ${B(o.L_a, 2)}`);
      R('Angle of twist', `φ = ${B(phi, 3)} rads`);
      R('Second derivative of φ', `φ'' = ${eng(o.phi2_deg / DEG)} rads/m^{2}`);
      R('Induced minor axis moment', `M_{yt} = M_{x} × φ / 1 rad = ${D(S.Myt, 1000, 2)} kNm`);
      R('Normal stress at flange tip due to M_{yt}', `σ_{byt} =M_{yt} / Z_{y} = ${D(S.sbyt, 1e6, 0)}${N2}`);
      R('Normal stress at flange tip due to warping', `σ_{w} = E_{S5950} × W_{n0} × φ'' / 1 rad = ${B(o.sig_w, 0)}${N2}`);
      R('Interaction index', `i_{b} = M_{x} × m_{LT} / M_{b} + (σ_{byt} + σ_{w}) / p_{y} × ( 1 + 0.5 × M_{x} × m_{LT} / M_{b}) = ${B(o.ib, 2)}`);
    }
    C(ok(x.buckling_ib) ? 'Pass - Combined bending and torsion check satisfied' : 'FAIL - Combined bending and torsion check not satisfied');
    H('Local capacity under combined bending & torsion');
    T('For simplicity, a conservative check is applied using the maximum stresses due to each of the separate load effects, even though these do not necessarily all occur at the same section along the member.');
    R('Max. direct stress due to M_{x}', `σ_{bx} = M_{x} / Z_{x} = ${D(S.sbx, 1e6, 0)}${N2}`);
    R('Combined stress - eqn 2.22', `${hr ? 'σ_{bx} + σ_{byt}' : 'σ_{bx} + σ_{byt} + σ_{w}'} = ${D(hr ? S.sbx + S.sbyt : S.sbx + S.sbyt + o.sig_w * 1e6, 1e6, 0)}${N2}`);
    R('Design strength', `p_{y} = ${B(o.py, 0)}${N2}`);
    C(ok(x.local) ? 'Pass - Local capacity' : 'FAIL - Local capacity');
    if (hr) {
      H('Combined shear stresses SCI-P-057 section 2.3');
      T('For simplicity, a conservative check is applied using the maximum shear stresses due to each of the separate load effects, even though these do not necessarily all occur at the same section along the member.');
      R('Max. shear stress due to bending', `τ_{bw} = F_{vy} × Q_{w} / (I_{x} × 2 × t) = ${B(o.tau_bw, 0)} N/mm^(2)`);
      R('Max. shear stresses due to torsion', `τ_{t} = T_{o} / C = ${B(o.tau_t, 0)} N/mm^(2)`);
      R('Amplified shear stress due to torsion', `τ_{vt} = τ_{t} × (1 + 0.5 × M_{x} × m_{LT} / M_{b}) = ${B(o.tau_vt, 0)} N/mm^(2)`);
      R('Combined shear due to bending & torsion', `τ = τ_{bw} + τ_{vt} = ${B(o.tau, 0)} N/mm^(2)`);
      R('Shear strength', `p_{v} = 0.6 × p_{y} = ${B(o.pv, 0)} N/mm^(2)`);
    } else {
      H('Combined shear stresses - SCI-P-057 section 2.3');
      T('For simplicity, a conservative check is applied using the maximum shear stresses due to each of the separate load effects, even though these do not necessarily all occur at the same section along the member.');
      R('Max shear stresses due to bending in web', `τ_{bw} = F_{vy} × Q_{w} / ( I_{x} × t ) = ${B(o.tau_bw, 0)}${N2}`);
      R('Max shear stresses due to bending in flange', `τ_{bf} = F_{vy} × Q_{f} / (I_{x} × T) = ${B(o.tau_bf, 0)}${N2}`);
      R('Max shear stresses due to torsion in web', `τ_{tw} = abs(G × t × φ' / 1rad) = ${B(o.tau_tw, 0)}${N2}`);
      R('Max shear stresses due to torsion in flange', `τ_{tf} = abs(G × T × φ' / 1 rad) = ${B(o.tau_tf, 0)}${N2}`);
      R('Max shear stresses due to warping in flange', `τ_{wf} = abs( -E_{S5950} × S_{w1} × φ''' / 1 rad / T) = ${B(o.tau_wf, 0)}${N2}`);
      R('Amp shear stress torsion & warping in web', `τ_{vtw} = τ_{tw} × (1 + 0.5 × M_{x} × m_{LT} / M_{b}) = ${B(o.tau_vtw, 0)}${N2}`);
      R('Amp shear stress torsion & warping in flange', `τ_{vtf} = (τ_{tf} + τ_{wf}) × (1 + 0.5 × M_{x} × m_{LT} / M_{b}) = ${B(o.tau_vtf, 0)}${N2}`);
      H('Combined shear stresses due to bending, torsion & warping:');
      R('Combined shear stresses in web', `τ_{w} = τ_{bw} + τ_{vtw} = ${B(o.tau_w, 0)}${N2}`);
      R('Combined shear stresses in flange', `τ_{f} = τ_{bf} + τ_{vtf} = ${B(o.tau_f, 0)}${N2}`);
      R('Shear strength', `p_{v} = 0.6 × p_{y} = ${B(o.pv, 0)}${N2}`);
    }
    C(ok(x.shear_stress) ? 'Pass - Combined shear stresses' : 'FAIL - Combined shear stresses');   // FAIL: live reference runs V1, V2
    if ('twist' in x) {
      H('Twist check');
      R('Total applied torque (unfactored)', `T_{qu} = ${D(S.Tqu, 1000, 2)} kNm`);
      if (p.twist_check === 'Lever') {
        R('Maximum twist under sls loading', `φ_{sls} = φ × T_{qu} / T_{q} = ${B(o.phi_sls_deg, 2)} degs`);
        R('Lever arm for deflection due to twist', `h_{δ} = ${B(p.h_d, 0)} mm`);
        R('Deflection due to twist', `δ_{h.sls} = h_{δ} × φ_{sls} / 1 rad = ${B(o.delta_h_sls, 1)} mm`);
        R('Deflection limit', `δ_{h.lim} = ${B(p.d_h_lim, 0)} mm`);
        C(ok(x.twist) ? 'Pass - Deflection due to twist' : 'FAIL - Deflection due to twist');   // live reference runs V2, V5
      } else {
        R('Maximum twist under sls loading', `φ_{sls} = φ × T_{qu} / T_{q} = ${B(o.phi_sls_deg, 2)} deg`);
        R('Twist limit', `φ_{lim} = ${B(p.phi_lim, 2)} deg`);
        C(ok(x.twist) ? 'Pass - Twist' : 'FAIL - Twist');
      }
    }
    if ('deflection' in x) {
      H('Deflection');
      R('Maximum y-axis deflection', `δ_{y_max} = ${D(S.dy, 1e-3, 1)} mm`);
      R('Deflection limit - cl. 2.5.2', `δ_{lim} =  min(L/ k_{δ}, δ_{lim_abs}) = ${D(S.dlimb, 1e-3, 1)} mm`);
      C(ok(x.deflection) ? 'Pass - Deflection within specified limit' : 'FAIL - Deflection exceeds specified limit');
    }
    return items;
  }

  // plain-text lines exactly as the reference copies them (for tests): rows "label;<TAB>expr", second tab stop doubled
  function lines(items) {
    const out = [];
    for (const it of items) {
      if (it.t === 'sketch' || it.t === 'ver') continue;
      if (it.t === 'row') out.push(['row', plain(it.l), plain(it.e)]);
      else if (it.t === 'row2') out.push(['row2', plain(it.l), plain(it.e)]);
      else out.push(['text', '', plain(it.s)]);
    }
    return out;
  }

  const SMSReport = { build, lines, plain, fx, sci, eng, RESTRAINT };
  if (typeof module !== 'undefined' && module.exports) module.exports = SMSReport;
  else root.SMSReport = SMSReport;
})(typeof window !== 'undefined' ? window : globalThis);

/* Cross-section sketch for the masonry support report (SVG fragment, pt units).
   Conventions: outer leaf / slab edge on the left, torsion beam (or wall) on the
   right, steel light grey (#D3D3D3) with hairline outlines, masonry and slabs stippled grey, load arrows at the
   load positions, dimensions with filled arrowheads and 6 pt figures: leaf widths, cavity and e_mb along the
   top, b_m and e_1 along the bottom, fixing dimensions on the right.
   sketch(p, o) -> SVG markup for a box 536 x 250 pt (x from 46, y from 0), centred on the body.
   gap(p)       -> distance (pt) from the version line to the heading after the sketch.
   the reference sizes its picture from its own drawing (RTF: scaled to floor(5040 twips / natural height) %, heading at
   0.6 x (round(height / 12 twips) + 56) pt): 280.8-285.6 pt in the 14 live documents. Tall sketches lose more to
   the floor - 280.8 / 281.4 for the channel and UKB 305 examples, 283.8-285.6 for the others - so the gap is taken
   from the sketch's own height/width: 281.4 pt from 1.17 up, else 284.4 (within 1.2 pt of every the reference print). */
(function (root) {
  'use strict';
  // placement measured on the 14 live the reference prints: the figure (dimension figures included) is 243.9 +- 1.2 pt tall,
  // starts 6.2 pt below the version line (0.2 pt into this fragment, which the renderer puts 6 pt below it) and is
  // centred at x = 300.2 +- 0.8 pt; widest figure 277.7 pt
  // dimensioning as the reference: everything scales with its canvas margin M = 0.14 x min(drawing width, height) - dimension
  // lines 0.47 M off the drawing, ticks and arrowheads 0.22 M, figures 0.30 M high
  const CX = 300.2, HT = 243.9, WMAX = 280, TOP = 0.2, K_OFF = 0.47, K_TICK = 0.22, K_TXT = 0.30;
  const STEEL = 'fill="#D3D3D3" stroke="#000" stroke-width="0.3" stroke-linejoin="round"';
  const MAS = 'fill="url(#cs-mas)" stroke="#8c8c8c" stroke-width="0.2"';
  const f0 = x => String(Math.round(x));

  function gap(p) {
    const { tor, xmin, xmax, ymin, ymax } = frame(p);
    return tor && (ymax - ymin) / (xmax - xmin) >= 1.17 ? 281.4 : 284.4;   // fixed to a wall: 283.8 / 285.6 in the reference
  }

  // model geometry (mm) and extents shared by sketch() and gap()
  function frame(p) {
    const tor = p.design_method === 'Masonry support and torsion', m = p.sb_method, bolt = p.fixing === 'Bolt';
    const slabS = p.sb_load === 'Slab' && m !== 'SCI', slabM = tor && p.mb_load === 'Slab';
    const s = bolt ? (p.shim || 0) : 0, emb = tor ? (p.e_mb || 0) : 0;
    let D = 0, B = 0, mb = p.mb;
    if (tor) { D = mb.D; B = mb.B; }
    // support member geometry (mm): x outward from the beam face (or wall face), y up from the beam soffit
    let lh, lv = 0, t, r = 0, ybear;
    if (m === 'Plate') { lh = p.lh; t = p.tsb; ybear = 0; }
    else if (m === 'Cantilever') { const lg = (p.orientation || 'Long') === 'Long'; lh = lg ? p.angle.D : p.angle.B; lv = lg ? p.angle.B : p.angle.D; t = p.angle.t; r = p.angle.r1 || 0; ybear = t; }
    else { lh = p.lh; lv = p.lv; t = p.tsb; r = p.r1sb; ybear = t; }
    const hb = tor && m === 'Cantilever' ? (p.h_bf || 0) : 0;          // angle soffit relative to the beam soffit (- below)
    const y0 = m === 'Plate' ? -t : hb;                                  // soffit of the support member
    const heel = m === 'Plate' ? 0 : s, toe = heel + lh;
    const xo = p.c - emb;                                               // inner face of the outer leaf
    const hIn = tor ? B / 2 : 0;                                         // inner leaf drawn B/2 high (the reference drawHeight = D + B/2)
    // fixed to a wall: the reference draws on a 304.8 mm canvas - masonry to 1.2 x 304.8, a slab wall to 304.8 + t (X04, X06)
    const topY = tor ? D + (slabM ? p.d_smb : hIn) : slabS ? 304.8 + t : 1.2 * 304.8;
    const outerTop = slabS ? (y0 + ybear + (m === 'Plate' ? t : 0)) + p.d_ssb : topY;
    // extents
    const wallW = slabS ? 60 : 0.3 * 304.8;                              // wall drawn 60 / 91.4 mm thick (X04 / X06)
    const xmin = tor ? -B - (slabM ? 40 : Math.max(0, emb + p.b_mmb - B)) : -s - wallW;
    const xmax = slabS ? Math.max(toe, heel + p.e1) + 60 : Math.max(toe, xo + p.b_msb);
    const ymin = tor ? Math.min(y0, 0) : Math.min(y0, 0) - 52, ymax = Math.max(topY, outerTop);   // wall runs 52 mm below
    return { tor, m, bolt, slabS, slabM, s, emb, D, B, mb, lh, lv, t, r, ybear, hb, y0, heel, toe, xo, hIn, topY, outerTop, xmin, xmax, ymin, ymax, wallW };
  }

  function sketch(p, o) {
    const { tor, m, bolt, slabS, slabM, s, emb, D, B, mb, lh, lv, t, r, ybear, hb, y0, heel, toe, xo, hIn, topY, outerTop, xmin, xmax, ymin, ymax, wallW } = frame(p);
    const M = 0.14 * Math.min(xmax - xmin, ymax - ymin);
    const sc = Math.min(HT / (ymax - ymin + 2 * (K_OFF + K_TXT * 0.62) * M), WMAX / (xmax - xmin + (K_OFF + K_TXT) * M));
    const FSD = K_TXT * M * sc, TICK = K_TICK * M * sc, AL = K_TICK * M * sc, AW = AL / 2;
    const midx = (xmin + xmax) / 2;
    const X = x => (CX - (x - midx) * sc), Y = y => TOP + (K_OFF + K_TXT * 0.62) * M * sc + (ymax - y) * sc;
    const P = (x, y) => `${X(x).toFixed(2)},${Y(y).toFixed(2)}`;
    const g = [];
    g.push('<defs><pattern id="cs-mas" patternUnits="userSpaceOnUse" width="1.6" height="1.6"><rect width="1.6" height="1.6" fill="#c8c8c8"/><circle cx="0.4" cy="0.4" r="0.32" fill="#9b9b9b"/><circle cx="1.2" cy="1.2" r="0.32" fill="#a8a8a8"/></pattern></defs>');
    const rect = (x1, y1, x2, y2, st) => g.push(`<polygon points="${P(x1, y1)} ${P(x2, y1)} ${P(x2, y2)} ${P(x1, y2)}" ${st}/>`);
    // ---- masonry / slab / wall
    if (!tor) rect(-s - wallW, ymin, -s, topY, MAS);
    if (slabS) rect(heel + t + 2, y0 + (m === 'Plate' ? t : t), Math.max(toe, heel + p.e1) + 60, outerTop, MAS);
    else rect(xo, y0 + (m === 'Plate' ? t : t), xo + p.b_msb, topY, MAS);
    if (tor) {
      if (slabM) rect(-B - 40, D, -emb, D + p.d_smb, MAS);
      else rect(-emb - p.b_mmb, D, -emb, D + hIn, MAS);
    }
    // ---- torsion beam
    if (tor) {
      const tw = mb.t, T = mb.T || mb.t;
      if (mb.type === 'C' && !p.mirror) g.push(`<polygon points="${P(0, 0)} ${P(-B, 0)} ${P(-B, T)} ${P(-tw, T)} ${P(-tw, D - T)} ${P(-B, D - T)} ${P(-B, D)} ${P(0, D)}" ${STEEL}/>`);
      else if (mb.type === 'C') g.push(`<polygon points="${P(-B, 0)} ${P(0, 0)} ${P(0, T)} ${P(-B + tw, T)} ${P(-B + tw, D - T)} ${P(0, D - T)} ${P(0, D)} ${P(-B, D)}" ${STEEL}/>`);
      else if (mb.type === 'I') {
        const a = -B / 2 - tw / 2, b = -B / 2 + tw / 2;
        g.push(`<polygon points="${P(-B, 0)} ${P(0, 0)} ${P(0, T)} ${P(b, T)} ${P(b, D - T)} ${P(0, D - T)} ${P(0, D)} ${P(-B, D)} ${P(-B, D - T)} ${P(a, D - T)} ${P(a, T)} ${P(-B, T)}" ${STEEL}/>`);
      } else {
        const ro = Math.min(2 * tw, B / 4) * sc, ri = Math.max(ro - tw * sc, 0.5);
        g.push(`<rect x="${X(0).toFixed(2)}" y="${Y(D).toFixed(2)}" width="${(B * sc).toFixed(2)}" height="${(D * sc).toFixed(2)}" rx="${ro.toFixed(2)}" ${STEEL}/>`);
        g.push(`<rect x="${X(-tw).toFixed(2)}" y="${Y(D - tw).toFixed(2)}" width="${((B - 2 * tw) * sc).toFixed(2)}" height="${((D - 2 * tw) * sc).toFixed(2)}" rx="${ri.toFixed(2)}" fill="#fff" stroke="#000" stroke-width="0.3"/>`);
      }
    }
    // ---- support member
    if (m === 'Plate') rect(-(p.lplate - p.lh), -t, p.lh, 0, STEEL);
    else {
      const bx = heel, by = y0;
      const R = r + t, pts = [P(bx, by + lv)];
      if (r > 0) {
        for (let i = 0; i <= 8; i++) { const a = Math.PI + i / 8 * Math.PI / 2; pts.push(P(bx + R + R * Math.cos(a), by + R + R * Math.sin(a))); }
        pts.push(P(bx + lh, by), P(bx + lh, by + t));
        for (let i = 0; i <= 8; i++) { const a = 1.5 * Math.PI - i / 8 * Math.PI / 2; pts.push(P(bx + R + r * Math.cos(a), by + R + r * Math.sin(a))); }
        pts.push(P(bx + t, by + lv));
      } else pts.push(P(bx, by), P(bx + lh, by), P(bx + lh, by + t), P(bx + t, by + t), P(bx + t, by + lv));
      g.push(`<polygon points="${pts.join(' ')}" ${STEEL}/>`);
      if (s > 0) rect(0, by, s, by + lv, 'fill="#fff" stroke="#000" stroke-width="0.3"');
      if (bolt) {
        const yb = by + lv - p.eb, d = Number(p.bolt_size.slice(1));
        rect(heel - 2, yb - d / 2, heel + t + 2, yb + d / 2, 'fill="#fff" stroke="#000" stroke-width="0.25"');
        g.push(`<line x1="${X(heel + t + 0.9 * d).toFixed(2)}" y1="${Y(yb).toFixed(2)}" x2="${X(tor ? -Math.min(mb.type === 'C' ? mb.t + 12 : B / 2, 40) : -45).toFixed(2)}" y2="${Y(yb).toFixed(2)}" stroke="#000" stroke-width="0.9"/>`);
      } else if (m === 'Cantilever' && tor) {
        const w = Math.max(p.s_weld || 6, 4);
        g.push(`<polygon points="${P(heel, by + lv)} ${P(heel + w, by + lv)} ${P(heel, by + lv + w)}" fill="#000"/>`);
      }
    }
    if (m === 'Plate' && tor) {   // fillet welds at the beam soffit
      const w = Math.max(p.s_weld || 6, 4);
      g.push(`<polygon points="${P(0, 0)} ${P(w, 0)} ${P(0, w)}" fill="#000"/>`, `<polygon points="${P(-B, 0)} ${P(-B - w, 0)} ${P(-B, w)}" fill="#000"/>`);
    }
    // ---- load arrows
    const arrow = (x, y1, y2) => {
      const xa = X(x), ya = Y(y1), yb = Y(y2);
      g.push(`<line x1="${xa.toFixed(2)}" y1="${ya.toFixed(2)}" x2="${xa.toFixed(2)}" y2="${(yb - 3).toFixed(2)}" stroke="#000" stroke-width="0.3"/>`,
        `<polygon points="${(xa - 1.6).toFixed(2)},${(yb - 3.2).toFixed(2)} ${(xa + 1.6).toFixed(2)},${(yb - 3.2).toFixed(2)} ${xa.toFixed(2)},${yb.toFixed(2)}" fill="#000"/>`);
    };
    const ySup = y0 + t;
    const xLoad = m === 'SCI' ? p.c - emb + (o.e || 0) : heel + p.e1;       // stainless: load e beyond the leaf face
    arrow(xLoad, ySup + Math.min(55, (outerTop - ySup) * 0.45), ySup);
    if (tor) arrow(slabM ? -p.e2 : -emb - p.b_mmb / 2, D + (slabM ? p.d_smb : hIn) * 0.75, D);
    // ---- dimensions
    const dimH = (x1, x2, y, label) => {
      if (Math.abs(x2 - x1) < 0.5) return;
      const a = Math.min(X(x1), X(x2)), b = Math.max(X(x1), X(x2)), yy = y;
      g.push(`<line x1="${a.toFixed(2)}" y1="${(yy - TICK).toFixed(2)}" x2="${a.toFixed(2)}" y2="${(yy + TICK).toFixed(2)}" stroke="#000" stroke-width="0.25"/>`,
        `<line x1="${b.toFixed(2)}" y1="${(yy - TICK).toFixed(2)}" x2="${b.toFixed(2)}" y2="${(yy + TICK).toFixed(2)}" stroke="#000" stroke-width="0.25"/>`);
      const tw = label.length * FSD * 0.556 + FSD / 2;
      const mid = (a + b) / 2, base = (yy + FSD * 0.36).toFixed(2);
      if (b - a <= tw + 2 * AL && b - a > label.length * FSD * 0.556 + 1) {      // narrow: arrows outside, figure inside
        g.push(`<polygon points="${a.toFixed(2)},${yy.toFixed(2)} ${(a - AL).toFixed(2)},${(yy - AW).toFixed(2)} ${(a - AL).toFixed(2)},${(yy + AW).toFixed(2)}" fill="#000"/>`,
          `<polygon points="${b.toFixed(2)},${yy.toFixed(2)} ${(b + AL).toFixed(2)},${(yy - AW).toFixed(2)} ${(b + AL).toFixed(2)},${(yy + AW).toFixed(2)}" fill="#000"/>`);
        g.push(`<text x="${mid.toFixed(2)}" y="${base}" font-size="${FSD.toFixed(2)}" text-anchor="middle">${label}</text>`);
        return;
      }
      if (b - a > tw + 2 * AL) {
        g.push(`<line x1="${a.toFixed(2)}" y1="${yy.toFixed(2)}" x2="${(mid - tw / 2).toFixed(2)}" y2="${yy.toFixed(2)}" stroke="#000" stroke-width="0.25"/>`,
          `<line x1="${(mid + tw / 2).toFixed(2)}" y1="${yy.toFixed(2)}" x2="${b.toFixed(2)}" y2="${yy.toFixed(2)}" stroke="#000" stroke-width="0.25"/>`);
        g.push(`<polygon points="${a.toFixed(2)},${yy.toFixed(2)} ${(a + AL).toFixed(2)},${(yy - AW).toFixed(2)} ${(a + AL).toFixed(2)},${(yy + AW).toFixed(2)}" fill="#000"/>`,
          `<polygon points="${b.toFixed(2)},${yy.toFixed(2)} ${(b - AL).toFixed(2)},${(yy - AW).toFixed(2)} ${(b - AL).toFixed(2)},${(yy + AW).toFixed(2)}" fill="#000"/>`);
        g.push(`<text x="${mid.toFixed(2)}" y="${base}" font-size="${FSD.toFixed(2)}" text-anchor="middle">${label}</text>`);
      } else g.push(`<text x="${mid.toFixed(2)}" y="${(yy - 2).toFixed(2)}" font-size="${FSD.toFixed(2)}" text-anchor="middle" transform="rotate(-90 ${mid.toFixed(2)} ${(yy - 2).toFixed(2)})">${label}</text>`);
    };
    const dimV = (y1, y2, x, label) => {
      if (Math.abs(y2 - y1) < 0.5) return;
      const a = Math.min(Y(y1), Y(y2)), b = Math.max(Y(y1), Y(y2));
      g.push(`<line x1="${(x - TICK).toFixed(2)}" y1="${a.toFixed(2)}" x2="${(x + TICK).toFixed(2)}" y2="${a.toFixed(2)}" stroke="#000" stroke-width="0.25"/>`,
        `<line x1="${(x - TICK).toFixed(2)}" y1="${b.toFixed(2)}" x2="${(x + TICK).toFixed(2)}" y2="${b.toFixed(2)}" stroke="#000" stroke-width="0.25"/>`);
      const th = label.length * FSD * 0.556 + FSD / 2, mid = (a + b) / 2;
      if (b - a > th + 2 * AL) {
        g.push(`<line x1="${x.toFixed(2)}" y1="${a.toFixed(2)}" x2="${x.toFixed(2)}" y2="${(mid - th / 2).toFixed(2)}" stroke="#000" stroke-width="0.25"/>`,
          `<line x1="${x.toFixed(2)}" y1="${(mid + th / 2).toFixed(2)}" x2="${x.toFixed(2)}" y2="${b.toFixed(2)}" stroke="#000" stroke-width="0.25"/>`,
          `<polygon points="${x.toFixed(2)},${a.toFixed(2)} ${(x - AW).toFixed(2)},${(a + AL).toFixed(2)} ${(x + AW).toFixed(2)},${(a + AL).toFixed(2)}" fill="#000"/>`,
          `<polygon points="${x.toFixed(2)},${b.toFixed(2)} ${(x - AW).toFixed(2)},${(b - AL).toFixed(2)} ${(x + AW).toFixed(2)},${(b - AL).toFixed(2)}" fill="#000"/>`);
      }
      g.push(`<text x="${(x + FSD * 0.36).toFixed(2)}" y="${mid.toFixed(2)}" font-size="${FSD.toFixed(2)}" text-anchor="middle" transform="rotate(-90 ${(x + FSD * 0.36).toFixed(2)} ${mid.toFixed(2)})">${label}</text>`);
    };
    const ytop = Y(ymax + K_OFF * M);
    if (!slabS) dimH(xo + p.b_msb, xo, ytop, f0(p.b_msb));
    dimH(xo, tor ? -emb : 0, ytop, f0(slabS ? p.c - emb : p.c - emb));
    if (tor && emb > 0) dimH(-emb, 0, ytop, f0(emb));
    if (tor && !slabM) dimH(-emb, -emb - p.b_mmb, ytop, f0(p.b_mmb));
    if (slabM) dimH(-emb, -p.e2, ytop, f0(p.e2 - emb));
    const ybot = Y(Math.min(y0, 0) - K_OFF * M);
    dimH(heel + lh, xLoad, ybot, f0(heel + lh - xLoad));
    dimH(xLoad, heel, ybot, f0(xLoad - heel));
    if (!tor) dimV(y0, y0 + t, X(Math.max(toe, xo + p.b_msb)) + 10, f0(t));          // the reference: support thickness when fixed to a wall
    if (s > 0) dimH(s, 0, ybot, f0(s));
    const xr = X(tor ? -B : -s - wallW) + K_OFF * M * sc;
    if (bolt) { dimV(y0 + lv, y0 + lv - p.eb, xr, f0(p.eb)); dimV(y0 + lv - p.eb, y0, xr, f0(o.lb || lv - p.eb)); if (hb) dimV(Math.min(0, y0), Math.max(0, y0), xr, f0(Math.abs(hb))); }
    else if (m === 'Cantilever') { dimV(y0 + lv, y0, xr, f0(lv)); if (hb) dimV(Math.min(0, y0), Math.max(0, y0), xr, f0(Math.abs(hb))); }
    else if (m === 'Plate') dimV(0, -t, X(p.lh) - 8, f0(t));
    if (slabS) dimV(ySup, outerTop, X(Math.max(toe, heel + p.e1) + 60) - 8, f0(p.d_ssb));
    if (slabM) dimV(D, D + p.d_smb, X(-B - 40) + 10, f0(p.d_smb));
    // centre the whole figure (dimension lines and figures included) on x = CX, as the reference centres its picture
    const svg = g.join(''), xs = [];
    for (const m of svg.matchAll(/points="([^"]+)"/g)) for (const pt of m[1].trim().split(/\s+/)) xs.push(+pt.split(',')[0]);
    for (const m of svg.matchAll(/ x[12]?="([-\d.]+)"/g)) xs.push(+m[1]);
    const dx = xs.length ? CX - (Math.min(...xs) + Math.max(...xs)) / 2 : 0;
    return `<g transform="translate(${dx.toFixed(2)} 0)">${svg}</g>`;
  }
  root.SheetSketch = { sketch, gap, frame };
})(typeof window !== 'undefined' ? window : globalThis);

/* Calculation-sheet pages: lays out the items of SMSReport.build() in the reference layout and returns
   one SVG per page (coordinates in pt, measured from the reference prints):
     header grid 46.8-582.0 x 22.9-106.3 pt, body frame 46.0-582.0 x 113.4-753.5 pt (Letter; A4 keeps the
     same layout and extends the frame), labels at x 68.88, expressions at 277.8 (second tab 303.0), check
     lines right-aligned at 577.9, Arial 8.98 pt, sub/superscripts 5.75 pt (-0.48 / +3.0 pt), header labels 6.95 pt.
   Line pitch (exact on all six worked examples): 13.68 pt + D(previous) + A(next), with
     plain line 0/0, line with a bold value +0.36/-0.36, line with Greek (Symbol font) +0.96/+0.36,
     heading +4.32/-0.36, check line +0.12/-0.12, title +9.36/-0.36 (A/D in pt).
   A page holds lines whose baseline + D <= 745.5 pt; headings keep with the next line; the first line of a
   page sits at 127.56 + A. Browser only (text widths are measured on a canvas). */
(function (root) {
  'use strict';
  const LEFT = 68.88, EXPR = 277.8, EXPR2 = 303.0, RIGHT = 577.9, WRAP = 578.0;
  const FS = 8.98, FSS = 5.75, FSH = 6.95, FSD = 5.98;
  const TOP = 127.56, PITCH = 13.68;
  const SYM_CHARS = 'αβγδεηθλμπρστφχψωΔ×√≤≥';
  const PUA = { 'α': 'a', 'β': 'b', 'γ': 'g', 'δ': 'd', 'ε': 'e', 'η': 'h', 'θ': 'q', 'λ': 'l', 'μ': 'm', 'π': 'p', 'ρ': 'r',
    'σ': 's', 'τ': 't', 'φ': 'f', 'χ': 'c', 'ψ': 'y', 'ω': 'w', 'Δ': 'D', '×': '´', '√': 'Ö', '≤': '£', '≥': '³' };
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  // Greek letters in the Symbol font where it is installed; elsewhere Unicode Greek in Times New Roman
  const useSymbolFont = typeof navigator !== 'undefined' && /Win/.test(navigator.platform || navigator.userAgent || '');
  const symGlyph = ch => useSymbolFont ? String.fromCharCode(0xF000 + PUA[ch].charCodeAt(0)) : ch;
  const SYMFAMILY = useSymbolFont ? 'Symbol' : "'Times New Roman', serif";

  // ADV-BEGIN (generated by advances.py --json + put_advances.py: mean advance in em of each 8.98 pt
  // glyph in the 14 live the reference prints, per font class r / b / bi / sym)
  const ADV = {
    "b5|-": 0.333007, "bi| ": 0.278306, "bi|-": 0.333011, "bi|.": 0.282061, "bi|1": 0.556153, "bi|A": 0.722167,
    "bi|B": 0.722167, "bi|C": 0.722166, "bi|D": 0.722167, "bi|F": 0.610841, "bi|I": 0.27783, "bi|L": 0.619445,
    "bi|M": 0.819999, "bi|P": 0.666992, "bi|S": 0.666992, "bi|T": 0.619998, "bi|a": 0.558538, "bi|b": 0.613511,
    "bi|c": 0.558265, "bi|d": 0.614321, "bi|e": 0.557509, "bi|f": 0.336435, "bi|g": 0.613061, "bi|h": 0.615143,
    "bi|i": 0.279596, "bi|k": 0.563536, "bi|l": 0.280455, "bi|m": 0.896666, "bi|n": 0.61606, "bi|o": 0.614753,
    "bi|p": 0.617963, "bi|q": 0.613127, "bi|r": 0.386465, "bi|s": 0.561018, "bi|t": 0.336991, "bi|u": 0.614766,
    "bi|v": 0.570002, "bi|w": 0.771804, "bi|x": 0.567555, "bi|y": 0.556152, "b| ": 0.278995, "b|&": 0.722168,
    "b|(": 0.335728, "b|)": 0.339999, "b|,": 0.277834, "b|-": 0.336139, "b|.": 0.280587, "b|0": 0.561854,
    "b|1": 0.558231, "b|2": 0.559488, "b|3": 0.560227, "b|4": 0.560587, "b|5": 0.561856, "b|6": 0.559949,
    "b|7": 0.561532, "b|8": 0.561342, "b|9": 0.563158, "b|:": 0.319999, "b|A": 0.710001, "b|B": 0.722169,
    "b|C": 0.722168, "b|D": 0.722168, "b|E": 0.665356, "b|F": 0.61084, "b|G": 0.770001, "b|I": 0.279445, "b|L": 0.616221,
    "b|M": 0.84, "b|N": 0.722168, "b|O": 0.775639, "b|P": 0.664221, "b|R": 0.722169, "b|S": 0.665727, "b|T": 0.612855,
    "b|U": 0.722168, "b|W": 0.943848, "b|Y": 0.666993, "b|a": 0.557483, "b|b": 0.615258, "b|c": 0.558023,
    "b|d": 0.615753, "b|e": 0.557352, "b|f": 0.337669, "b|g": 0.614175, "b|h": 0.615334, "b|i": 0.27972, "b|k": 0.565384,
    "b|l": 0.282763, "b|m": 0.894716, "b|n": 0.614787, "b|o": 0.614774, "b|p": 0.615121, "b|q": 0.61084, "b|r": 0.386991,
    "b|s": 0.560339, "b|t": 0.335896, "b|u": 0.614753, "b|v": 0.547857, "b|w": 0.800984, "b|x": 0.556153,
    "b|y": 0.507285, "r5|,": 0.270636, "r5|.": 0.277831, "r5|0": 0.556151, "r5|1": 0.556152, "r5|2": 0.556149,
    "r5|3": 0.556149, "r5|5": 0.556151, "r5|9": 0.569999, "r5|E": 0.666992, "r5|G": 0.770001, "r5|L": 0.561007,
    "r5|Q": 0.770002, "r5|S": 0.666992, "r5|U": 0.740002, "r5|_": 0.556151, "r5|a": 0.563896, "r5|b": 0.556877,
    "r5|c": 0.500004, "r5|d": 0.562131, "r5|e": 0.560749, "r5|f": 0.272151, "r5|h": 0.556156, "r5|i": 0.209996,
    "r5|l": 0.21233, "r5|m": 0.835838, "r5|n": 0.556156, "r5|o": 0.568543, "r5|p": 0.556151, "r5|q": 0.55615,
    "r5|s": 0.500501, "r5|t": 0.270001, "r5|v": 0.520003, "r5|w": 0.689999, "r5|x": 0.500001, "r5|y": 0.481727,
    "r| ": 0.281051, "r|&": 0.666992, "r|'": 0.18728, "r|(": 0.336468, "r|)": 0.337164, "r|+": 0.583984, "r|,": 0.278158,
    "r|-": 0.336965, "r|.": 0.280841, "r|/": 0.278699, "r|0": 0.557874, "r|1": 0.557519, "r|2": 0.56057, "r|3": 0.559487,
    "r|4": 0.56191, "r|5": 0.558373, "r|6": 0.558112, "r|7": 0.563542, "r|8": 0.561255, "r|9": 0.563318, "r|;": 0.277,
    "r|=": 0.583996, "r|A": 0.666992, "r|B": 0.668709, "r|C": 0.722604, "r|D": 0.722168, "r|E": 0.666993, "r|F": 0.61084,
    "r|G": 0.778256, "r|H": 0.722168, "r|I": 0.277832, "r|J": 0.510001, "r|K": 0.666989, "r|L": 0.556815,
    "r|M": 0.810149, "r|N": 0.722168, "r|O": 0.77, "r|P": 0.666026, "r|R": 0.722169, "r|S": 0.666648, "r|T": 0.601834,
    "r|U": 0.722169, "r|V": 0.659999, "r|W": 1.0, "r|X": 0.659998, "r|[": 0.277829, "r|]": 0.280873, "r|^": 0.450001,
    "r|a": 0.559279, "r|b": 0.561629, "r|c": 0.508198, "r|d": 0.559282, "r|e": 0.559836, "r|f": 0.279569, "r|g": 0.56179,
    "r|h": 0.559701, "r|i": 0.227922, "r|k": 0.509077, "r|l": 0.22866, "r|m": 0.843255, "r|n": 0.560546, "r|o": 0.55824,
    "r|p": 0.560206, "r|q": 0.557962, "r|r": 0.337339, "r|s": 0.509096, "r|t": 0.280233, "r|u": 0.558362,
    "r|v": 0.495899, "r|w": 0.707251, "r|x": 0.481747, "r|y": 0.498292, "r|z": 0.494102, "r|\u2013": 0.563076,
    "sym|\uf020": 0.25381, "sym|\uf064": 0.494269, "sym|\uf065": 0.440836, "sym|\uf066": 0.494269, "sym|\uf068": 0.60114,
    "sym|\uf06c": 0.561062, "sym|\uf070": 0.534343, "sym|\uf073": 0.587782, "sym|\uf074": 0.427475,
    "sym|\uf0b4": 0.548186, "sym|\uf0d6": 0.547702
  };
  // ADV-END
  let ctx = null;
  // Text width as the reference sets it: its glyph advances differ from the font's (8.98 pt: 0-4 % wider, 0.74 % on average;
  // 5.75 pt: -4 to +4 %), so body text and sub/superscripts are measured with the reference own mean advances (ADV) and
  // anything else with the font advances x 1.0074.
  const ADV_SCALE = 1.0074;
  function measure(text, size, bold, italic, sym) {
    if (!ctx) ctx = document.createElement('canvas').getContext('2d');
    ctx.font = `${italic ? 'italic ' : ''}${bold ? 'bold ' : ''}${size}pt ${sym ? SYMFAMILY : 'Arial'}`;
    const font = s => ctx.measureText(s).width * 0.75 * ADV_SCALE;     // px -> pt
    const small = Math.abs(size - FSS) < 0.01;                          // sub/superscripts: the reference 5.75 pt advances
    if (Math.abs(size - FS) > 0.01 && !small) return font(text);
    const cls = (sym ? 'sym' : ((bold ? 'b' : '') + (italic ? 'i' : '')) || 'r') + (small ? '5' : '');
    let w = 0;
    for (const ch of text) { const a = ADV[cls + '|' + ch]; w += a !== undefined ? a * size : font(ch); }
    return w;
  }

  // markup -> runs {t, b, sub, sup, sym, i}
  function runs(m, base) {
    base = base || {};
    const out = [];
    let bold = !!base.b, i = 0, buf = '';
    const flush = (extra) => { if (buf) { out.push(Object.assign({ t: buf, b: bold, i: !!base.i }, extra || {})); buf = ''; } };
    while (i < m.length) {
      if (m.startsWith('**', i)) { flush(); bold = !bold; i += 2; continue; }
      if ((m[i] === '_' || m[i] === '^') && m[i + 1] === '{') {
        flush();
        const j = m.indexOf('}', i + 2), inner = m.slice(i + 2, j), kind = m[i] === '_' ? 'sub' : 'sup';
        for (const ch of inner) {
          if (SYM_CHARS.includes(ch)) out.push({ t: ch, b: bold, sym: true, [kind]: true });
          else if (out.length && out[out.length - 1][kind] && !out[out.length - 1].sym && out[out.length - 1].b === bold) out[out.length - 1].t += ch;
          else out.push({ t: ch, b: bold, [kind]: true });
        }
        i = j + 1; continue;
      }
      const ch = m[i];
      if (SYM_CHARS.includes(ch)) { flush(); out.push({ t: ch, b: bold, sym: true, i: !!base.i }); i++; continue; }
      buf += ch; i++;
    }
    flush();
    return out;
  }
  const runWidth = (r, size) => measure(r.sym ? symGlyph(r.t) : r.t, r.sub || r.sup ? FSS : (size || FS), r.b, r.i, r.sym);
  const width = (rs, size) => rs.reduce((a, r) => a + runWidth(r, size), 0);

  // split runs into words (a word keeps its trailing spaces), for line breaking
  function words(rs) {
    const out = [];
    let cur = [];
    for (const r of rs) {
      if (r.sym || r.sub || r.sup) { cur.push(r); continue; }
      const parts = r.t.split(/(?<= )/);
      parts.forEach((part, k) => {
        cur.push(Object.assign({}, r, { t: part }));
        if (part.endsWith(' ') && (k < parts.length - 1 || true)) { out.push(cur); cur = []; }
      });
    }
    if (cur.length) out.push(cur);
    return out;
  }
  function wrap(rs, x0, xc) {
    const lines = [];
    let line = [], x = x0, start = x0;
    for (const w of words(rs)) {
      const ww = width(w), wtrim = width(w.map((r, i) => i === w.length - 1 ? Object.assign({}, r, { t: r.t.replace(/ +$/, '') }) : r));
      if (line.length && x + wtrim > WRAP) { lines.push({ x: start, runs: line }); line = []; x = xc; start = xc; }
      line.push(...w); x += ww;
    }
    if (line.length || !lines.length) lines.push({ x: start, runs: line });
    return lines;
  }
  // only full-size Greek raises the line (a Greek subscript such as h_delta keeps the plain/bold metrics: live run V2)
  const hasSym = rs => rs.some(r => r.sym && !r.sub && !r.sup), hasBold = rs => rs.some(r => r.b && r.t.trim());

  // ------------------------------------------------------------------ items -> lines
  function toLines(items) {
    const L = [];
    for (const it of items) {
      if (it.t === 'title') { L.push({ k: 'title', a: 9.36, d: -0.36, segs: [{ x: LEFT, runs: runs(it.s, { b: 1 }) }], keep: true, underline: true }); continue; }
      if (it.t === 'sub') { L.push({ k: 'head', a: 4.32, d: -0.36, segs: [{ x: LEFT, runs: runs(it.s, { b: 1 }) }], keep: true }); continue; }
      if (it.t === 'ver') { L.push({ k: 'ver', a: -1.68, d: -1.92, segs: [{ x: RIGHT, anchor: 'end', size: FSH, runs: [{ t: it.s }] }], keep: true }); continue; }
      if (it.t === 'sketch') { L.push({ k: 'sketch', keep: true }); continue; }
      if (it.t === 'h') {
        wrap(runs(it.s, { b: 1 }), LEFT, LEFT).forEach((seg, i, arr) => L.push({ k: 'head', a: i ? 0.36 : 4.32, d: -0.36, segs: [seg], keep: true, cont: i > 0 }));
        continue;
      }
      if (it.t === 'txt') {
        wrap(runs(it.s), LEFT, LEFT).forEach((seg, i) => L.push({ k: 'row', segs: [seg], cont: i > 0, glue: i > 0 }));
        continue;
      }
      if (it.t === 'chk') {
        const rs = runs(it.s, { b: 1, i: 1 });
        L.push({ k: 'chk', a: 0.12, d: -0.12, segs: [{ x: RIGHT, anchor: 'end', runs: rs }] });
        continue;
      }
      // rows
      const lab = it.l ? runs(it.l) : [];
      const semi = { t: ';', hidden: true };
      let x0 = it.t === 'row2' ? EXPR2 : EXPR;
      const ex = wrap(runs(it.e), x0, EXPR);
      ex.forEach((seg, i) => {
        const segs = [];
        if (i === 0 && it.l) segs.push({ x: LEFT, runs: lab.concat([semi]) });
        if (i === 0 && !it.l && it.t !== 'row2') segs.push({ x: LEFT, runs: [] });
        segs.push(seg);
        const all = segs.flatMap(s => s.runs).filter(r => !r.hidden);
        const g = hasSym(all), b = hasBold(all);
        L.push({ k: 'row', a: g ? 0.96 : b ? 0.36 : 0, d: g ? 0.36 : b ? -0.36 : 0, segs, glue: i > 0 });
      });
    }
    for (const l of L) if (l.k === 'row' && l.a === undefined) { const all = l.segs.flatMap(s => s.runs); const g = hasSym(all), b = hasBold(all); l.a = g ? 0.96 : b ? 0.36 : 0; l.d = g ? 0.36 : b ? -0.36 : 0; }
    return L;
  }

  // ------------------------------------------------------------------ pagination
  function paginate(L, bottom, sketchGap) {
    const pages = [[]];
    let y = null, prev = null;
    const fits = (yy, l) => yy + (l.d || 0) <= bottom;
    for (let i = 0; i < L.length; i++) {
      const l = L[i];
      if (l.k === 'sketch') { l.y = y; pages[pages.length - 1].push(l); prev = { d: 0, sketch: true, y: y }; continue; }
      let ny;
      if (y === null) ny = TOP + l.a;
      else if (prev && prev.sketch) ny = prev.y + (sketchGap || 284.4) - (4.32 - l.a);   // heading after the sketch: version line + gap (SheetSketch.gap)
      else ny = y + PITCH + prev.d + l.a;
      // keep headings / titles with the next line, and wrapped lines with their first line
      let need = ny, j = i, yy = ny, pl = l;
      while (L[j] && L[j].keep && L[j + 1] && L[j + 1].k !== 'sketch') { yy = yy + PITCH + pl.d + L[j + 1].a; pl = L[j + 1]; j++; need = yy; }
      while (L[j + 1] && L[j + 1].glue) { yy = yy + PITCH + L[j].d + L[j + 1].a; j++; need = yy; pl = L[j]; }
      const endFits = need + (pl.d || 0) <= bottom;
      if (y !== null && (!fits(ny, l) || !endFits) && !(prev && prev.sketch) && !l.glue) {
        pages.push([]);
        ny = TOP + l.a;
      }
      l.y = ny; y = ny; prev = l;
      pages[pages.length - 1].push(l);
    }
    return pages;
  }

  // ------------------------------------------------------------------ SVG
  function tspan(r, size, x) {
    const sz = r.sub || r.sup ? FSS : (size || FS);
    const attrs = [];
    if (x !== undefined) attrs.push(`x="${x.toFixed(2)}"`);
    if (r.hidden) attrs.push('fill="#fff"');
    if (r.b) attrs.push('font-weight="bold"');
    if (r.i) attrs.push('font-style="italic"');
    if (r.sym) attrs.push(`font-family="${SYMFAMILY}"`);
    if (sz !== (size || FS)) attrs.push(`font-size="${sz}"`);
    if (r.sub) attrs.push('baseline-shift="-0.48"');
    if (r.sup) attrs.push('baseline-shift="3"');
    return `<tspan ${attrs.join(' ')}>${esc(r.sym ? symGlyph(r.t) : r.t)}</tspan>`;
  }
  // every word (and every Greek letter, sub- or superscript) starts where the reference advances put it, so the glyphs line up
  // with the reference print along the whole line instead of drifting with the font's own advances
  function runsSvgAt(rs, x0, size) {
    let x = x0;
    const out = [];
    for (const r of rs) {
      for (const part of (r.sym || r.sub || r.sup) ? [r.t] : r.t.split(/(?<= )/)) {
        if (!part) continue;
        const w = Object.assign({}, r, { t: part });
        out.push(tspan(w, size, x));
        x += runWidth(w, size);
      }
    }
    return out.join('');
  }
  function header(meta, pageNo, W) {
    const ln = (x1, y1, x2, y2, w) => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="#010101" stroke-width="${w}"/>`;
    const tx = (x, y, s, size, extra) => `<text x="${x}" y="${y}" font-size="${size || FSH}"${extra || ''}>${esc(s)}</text>`;
    const h = [];
    // outer box and grid (line widths as the reference: 0.48 outer / horizontal, 0.24 inner verticals)
    h.push(ln(46.8, 22.92, 582.24, 22.92, 0.48), ln(46.8, 106.32, 582.24, 106.32, 0.48), ln(47.04, 22.68, 47.04, 106.56, 0.48), ln(582.0, 22.68, 582.0, 106.56, 0.48));
    h.push(ln(210.12, 22.68, 210.12, 106.56, 0.24), ln(465.36, 22.68, 465.36, 106.56, 0.24));
    h.push(ln(210.12, 50.64, 582.24, 50.64, 0.48), ln(210.12, 78.6, 582.24, 78.6, 0.48));
    for (const x of [273.96, 337.8, 401.52, 522.0]) h.push(ln(x, 78.6, x, 106.56, 0.24));
    // company block (no third-party logo)
    if (meta.company) h.push(tx(128.58, 51.36, meta.company, FS, ' text-anchor="middle"'));
    if (meta.address) h.push(tx(128.58, 63.36, meta.address, FSH, ' text-anchor="middle"'));
    if (meta.logo) h.push(`<image href="${esc(meta.logo)}" x="80" y="24.5" width="97" height="20" preserveAspectRatio="xMidYMid meet"/>`);
    const L = [[212.88, 31.80, 'Project', meta.project], [468.11, 31.80, 'Job Ref.', meta.jobref], [212.88, 59.52, 'Section', meta.section],
      [468.06, 59.52, 'Sheet no./rev.', ' ' + pageNo + (meta.rev ? '/' + meta.rev : '')], [212.88, 87.48, 'Calc. by', meta.calcby],
      [276.72, 87.48, 'Date', meta.calcdate], [340.56, 87.48, "Chk'd by", meta.chkby], [404.24, 87.48, 'Date', meta.chkdate],
      [468.14, 87.48, "App'd by", meta.appby], [524.70, 87.48, 'Date', meta.appdate]];
    const cellR = { 212.88: 465.36, 468.11: 582.0, 468.06: 582.0, 276.72: 337.8, 340.56: 401.52, 404.24: 465.36, 468.14: 522.0, 524.70: 582.0 };
    L.forEach(([x, y, lab, val], i) => {
      h.push(tx(x, y, lab, FSH));
      if (!val) return;
      const right = x === 212.88 && y === 87.48 ? 273.96 : cellR[x], id = `cs-clip-${pageNo}-${i}`;
      h.push(`<clipPath id="${id}"><rect x="${(x - 2.5).toFixed(2)}" y="${(y - 9).toFixed(2)}" width="${(right - x + 2).toFixed(2)}" height="27"/></clipPath>`);
      h.push(`<g clip-path="url(#${id})">${tx(x, +(y + 12.36).toFixed(2), val, FS)}</g>`);
    });
    return h.join('');
  }
  function pageSvg(lines, meta, pageNo, paper, sketchSvg) {
    const W = paper.w, H = paper.h, frameBottom = paper.frameBottom;
    const out = [];
    out.push(`<svg xmlns="http://www.w3.org/2000/svg" class="tpage" viewBox="0 0 ${W} ${H}" width="${W}pt" height="${H}pt" font-family="Arial" font-size="${FS}" fill="#000">`);
    out.push(`<rect x="0" y="0" width="${W}" height="${H}" fill="#fff"/>`);
    out.push(header(meta, pageNo, W));
    out.push(`<rect x="45.96" y="113.4" width="${(582.0 - 45.96).toFixed(2)}" height="${(frameBottom - 113.4).toFixed(2)}" fill="none" stroke="#010101" stroke-width="0.48"/>`);
    for (const l of lines) {
      if (l.k === 'sketch') { if (sketchSvg) out.push(`<g transform="translate(0 ${(l.y + 6).toFixed(2)})">${sketchSvg}</g>`); continue; }
      for (const s of l.segs) {
        if (!s.runs.length) continue;
        const size = s.size ? ` font-size="${s.size}"` : '';
        const x0 = s.anchor === 'end' ? s.x - width(s.runs, s.size) : s.x;     // right-aligned: start = right edge - the reference width
        out.push(`<text x="${x0.toFixed(2)}" y="${l.y.toFixed(2)}"${size} xml:space="preserve">${runsSvgAt(s.runs, x0, s.size)}</text>`);
        if (l.underline) out.push(`<line x1="${s.x.toFixed(2)}" y1="${(l.y + 2.38).toFixed(2)}" x2="${(s.x + width(s.runs)).toFixed(2)}" y2="${(l.y + 2.38).toFixed(2)}" stroke="#000" stroke-width="0.48"/>`);
      }
    }
    out.push('</svg>');
    return out.join('');
  }

  const PAPER = { Letter: { w: 612, h: 792, frameBottom: 753.48, bottom: 745.5 }, A4: { w: 595.28, h: 841.89, frameBottom: 803.37, bottom: 795.4 } };
  function render(items, meta, opts) {
    opts = opts || {};
    const paper = PAPER[opts.paper] || PAPER.Letter;
    const L = toLines(items);
    const pages = paginate(L, paper.bottom, opts.sketchGap);
    return pages.map((ls, i) => pageSvg(ls, meta || {}, i + 1, paper, opts.sketch || ''));
  }

  root.SheetRender = { render, toLines, paginate, runs, PAPER };
})(typeof window !== 'undefined' ? window : globalThis);

/* Steel masonry support (BS 5950) — small engines, JavaScript port of engine/sms (Python).
   Units N, mm; angles in radians inside. Works in the browser (window.SMS) and as a module (module.exports). */
(function (root) {
  'use strict';
  // The benchmark results were computed with pi = 3.1415926536 (lambda_LT implied by p_E and every degree
  // conversion sit 3.25e-12 below Math.PI values); the method uses it. Our own section geometry keeps Math.PI.
  const PI_T = 3.1415926536;
  const DEG = 180 / PI_T;

  // ---------------------------------------------------------------- materials
  const materials = (() => {
    const E = 205000, G = 78800, RHO = 78.5e-6;
    const PY = {
      S275: [[16, 275], [40, 265], [63, 255], [80, 245], [100, 235], [150, 225]],
      S355: [[16, 355], [40, 345], [63, 335], [80, 325], [100, 315], [150, 295]],
      S460: [[16, 460], [40, 440], [63, 430], [80, 410], [100, 400]],
    };
    const PW = { S275: { 35: 220, 42: 220, 50: 220 }, S355: { 35: 220, 42: 250, 50: 250 }, S460: { 35: 220, 42: 250, 50: 280 } };
    const AS = { M12: 84.3, M16: 157, M20: 245, M22: 303, M24: 353, M27: 459, M30: 561, M33: 694, M36: 817 };
    const PS = { '4.6': 160, '8.8': 375 }, PT = { '4.6': 195, '8.8': 450 };
    // SCI AD187 0.2 % proof stress: cold-rolled strip up to 6 mm, hot-rolled from 8 mm (as the reference; 7 mm taken as hot [verify])
    const Y02 = { '304': 230, '304L': 220, '316': 240, '316L': 240 }, Y02_HOT = { '304': 210, '304L': 200, '316': 220, '316L': 220 };
    const py = (grade, t) => { for (const [lim, v] of PY[grade]) if (t <= lim) return v; throw new Error(`${grade}: t=${t} outside Table 9`); };
    const pw = (grade, cls) => PW[grade][cls];
    const rnd10 = x => Math.floor(x / 10 + 0.5) * 10;
    const bolt = (grade, size) => ({ As: AS[size], Ps: rnd10(PS[grade] * AS[size]), Pt: rnd10(PT[grade] * AS[size]) });
    const stainless = (grade, gm = 1, t = 0) => (t <= 6 ? Y02 : Y02_HOT)[grade] / gm;
    return { E, G, RHO, py, pw, bolt, stainless };
  })();

  // ---------------------------------------------------------------- sections
  const sections = (() => {
    const K = (10 - 3 * Math.PI) / (12 - 3 * Math.PI);
    const SPA = 1 - Math.PI / 4;
    const SPI0 = 1 / 3 - (Math.PI / 16 - 4 / (9 * Math.PI)) - Math.PI / 4 * (1 - 4 / (3 * Math.PI)) ** 2;
    const SPIC = SPI0 - SPA * K * K;
    function iSection(D, B, t, T, r) {
      const a = SPA * r * r, ys = D / 2 - T - K * r, xs = t / 2 + K * r;
      const A = 2 * B * T + (D - 2 * T) * t + 4 * a;
      const Ix = 2 * (B * T ** 3 / 12 + B * T * (D / 2 - T / 2) ** 2) + t * (D - 2 * T) ** 3 / 12 + 4 * (SPIC * r ** 4 + a * ys * ys);
      const Iy = 2 * T * B ** 3 / 12 + (D - 2 * T) * t ** 3 / 12 + 4 * (SPIC * r ** 4 + a * xs * xs);
      const Sx = 2 * (B * T * (D - T) / 2 + t * (D / 2 - T) ** 2 / 2 + 2 * a * ys);
      const h = D - T;
      const a1 = -0.042 + 0.2204 * t / T + 0.1355 * r / T - 0.0865 * t * r / T ** 2 - 0.0725 * (t / T) ** 2;
      const D1 = ((T + r) ** 2 + (r + 0.25 * t) * t) / (2 * r + T);
      const J = 2 / 3 * B * T ** 3 + (D - 2 * T) * t ** 3 / 3 + 2 * a1 * D1 ** 4 - 0.42 * T ** 4;
      const H = Iy * h * h / 4, gam = 1 - Iy / Ix;
      return { type: 'I', D, B, t, T, r, A, Ix, Iy, Zx: Ix / (D / 2), Zy: Iy / (B / 2), Sx, ry: Math.sqrt(Iy / A), J, H,
        u: (4 * Sx * Sx * gam / (A * A * h * h)) ** 0.25, x: 0.566 * h * Math.sqrt(A / J), cy: 0 };
    }
    function channel(D, B, t, T, r, JH) {
      const a = SPA * r * r, ys = D / 2 - T - K * r;
      const parts = [[B * T, B / 2, D / 2 - T / 2, B * T ** 3 / 12, T * B ** 3 / 12], [B * T, B / 2, -(D / 2 - T / 2), B * T ** 3 / 12, T * B ** 3 / 12],
        [(D - 2 * T) * t, t / 2, 0, t * (D - 2 * T) ** 3 / 12, (D - 2 * T) * t ** 3 / 12],
        [a, t + K * r, ys, SPIC * r ** 4, SPIC * r ** 4], [a, t + K * r, -ys, SPIC * r ** 4, SPIC * r ** 4]];
      const A = parts.reduce((s, p) => s + p[0], 0), cy = parts.reduce((s, p) => s + p[0] * p[1], 0) / A;
      const Ix = parts.reduce((s, p) => s + p[3] + p[0] * p[2] * p[2], 0), Iy = parts.reduce((s, p) => s + p[4] + p[0] * (p[1] - cy) ** 2, 0);
      const Sx = 2 * (B * T * (D - T) / 2 + t * (D / 2 - T) ** 2 / 2 + a * ys), h = D - T;
      const a3 = -0.0908 + 0.2621 * t / T + 0.1231 * r / T - 0.0752 * t * r / T ** 2 - 0.0945 * (t / T) ** 2;
      const D3 = 2 * ((3 * r + t + T) - Math.sqrt(2 * (2 * r + t) * (2 * r + T)));
      let J = 2 / 3 * B * T ** 3 + (D - 2 * T) * t ** 3 / 3 + 2 * a3 * D3 ** 4 - 0.42 * T ** 4;
      let H = h * h / 4 * (Iy - A * (cy - t / 2) ** 2 * (h * h * A / (4 * Ix) - 1));
      if (JH) { J = JH[0]; H = JH[1]; }
      const gam = 1 - Iy / Ix;
      return { type: 'C', D, B, t, T, r, A, Ix, Iy, Zx: Ix / (D / 2), Zy: Iy / (B - cy), Sx, ry: Math.sqrt(Iy / A), J, H,
        u: (Iy * Sx * Sx * gam / (A * A * H)) ** 0.25, x: 1.132 * Math.sqrt(A * H / (Iy * J)), cy };
    }
    function rhsHot(D, B, t) {
      const solid = (d, b, rr) => { const a = SPA * rr * rr, yc = d / 2 - K * rr, xc = b / 2 - K * rr;
        return [b * d - 4 * a, b * d ** 3 / 12 - 4 * (SPIC * rr ** 4 + a * yc * yc), d * b ** 3 / 12 - 4 * (SPIC * rr ** 4 + a * xc * xc), b * d * d / 4 - 4 * a * yc]; };
      const o = solid(D, B, 1.5 * t), i = solid(D - 2 * t, B - 2 * t, t);
      const A = o[0] - i[0], Ix = o[1] - i[1], Iy = o[2] - i[2], Sx = o[3] - i[3];
      const Rc = 1.25 * t, hm = 2 * ((B - t) + (D - t)) - 2 * Rc * (4 - Math.PI), Ah = (B - t) * (D - t) - Rc * Rc * (4 - Math.PI);
      const Kc = 2 * Ah * t / hm, J = t ** 3 * hm / 3 + 2 * Kc * Ah;
      return { type: 'HR', D, B, t, T: t, r: 0, A, Ix, Iy, Zx: Ix / (D / 2), Zy: Iy / (B / 2), Sx, ry: Math.sqrt(Iy / A), J, C: J / (t + Kc / t) };
    }
    function angle(aL, bS, t, r1, r2) {
      const ar1 = SPA * r1 * r1, ar2 = SPA * r2 * r2;
      const A = aL * t + (bS - t) * t + ar1 - 2 * ar2;
      const Sy0 = aL * t * t / 2 + (bS - t) * t * (t + (bS - t) / 2) + ar1 * (t + K * r1) - ar2 * (t - K * r2) - ar2 * (bS - K * r2);
      const Sx0 = aL * t * aL / 2 + (bS - t) * t * t / 2 + ar1 * (t + K * r1) - ar2 * (aL - K * r2) - ar2 * (t - K * r2);
      return { A, cx: Sx0 / A, cy: Sy0 / A };
    }
    function torsionOpen(mb) {
      const { D, B, t, T, r } = mb, af = SPA * r * r, yf = D / 2 - T - K * r, h1 = D - T;
      if (mb.type === 'C') {
        const b1 = B - t / 2, e0 = 3 * b1 * b1 * T / (6 * b1 * T + h1 * t);
        return { e0, Wn0: h1 / 2 * (b1 - e0), Sw1: T * (b1 - e0) ** 2 * h1 / 4, Qf: (B - t) * T * (D - T) / 2 + af * yf, Qw: mb.Sx / 2 };
      }
      return { e0: 0, Wn0: h1 * B / 4, Sw1: h1 * B * B * T / 16, Qf: (B - t) / 2 * T * (D - T) / 2 + af * yf, Qw: mb.Sx / 2 };
    }
    function classify(mb, py) {
      const eps = Math.sqrt(275 / py);
      let bt, dt, fl, wb;
      if (mb.type === 'HR') { bt = (mb.B - 3 * mb.t) / mb.t; dt = (mb.D - 3 * mb.t) / mb.t; fl = [28, 32, 40]; wb = [64, 80, 120]; }
      else { bt = mb.B / 2 / mb.T; dt = (mb.D - 2 * mb.T - 2 * mb.r) / mb.t; fl = [9, 10, 15]; wb = [80, 100, 120]; }
      const c = (v, lims) => { const i = lims.findIndex(l => v <= l * eps); return i < 0 ? 3 : i; };
      return { cls: ['Plastic', 'Compact', 'Semi-Compact', 'Slender'][Math.max(c(bt, fl), c(dt, wb))], bt, dt };
    }
    // UKPFC torsion constants (J mm4, H mm6) as held in the reference 2019 data list
    const UKPFC_JH = {
      'UKPFC 100x50x10': [24953.00759, 496101000.0],
      'UKPFC 125x65x15': [46565.4754, 1979285000.0],
      'UKPFC 150x75x18': [60267.93592, 4729098000.0],
      'UKPFC 150x90x24': [116829.2342, 8964160000.0],
      'UKPFC 180x75x20': [72585.45186, 7620874000.0],
      'UKPFC 180x90x26': [131795.3553, 14217730000.0],
      'UKPFC 200x75x23': [109629.8507, 10823890000.0],
      'UKPFC 200x90x30': [181071.147, 19732160000.0],
      'UKPFC 230x75x26': [116841.1886, 15341900000.0],
      'UKPFC 230x90x32': [191148.1569, 27925290000.0],
      'UKPFC 260x75x28': [116052.4632, 20391470000.0],
      'UKPFC 260x90x35': [203745.9197, 37936610000.0],
      'UKPFC 300x100x46': [364496.8512, 81470240000.0],
      'UKPFC 300x90x41': [284972.4446, 57886930000.0],
      'UKPFC 380x100x54': [451929.4164, 150332100000.0],
      'UKPFC 430x100x64': [623771.3738, 217479100000.0]
    };
    function ownProps(mb) {
      if (mb.type === 'I') return Object.assign(iSection(mb.D, mb.B, mb.t, mb.T, mb.r), { name: mb.name });
      if (mb.type === 'C') return Object.assign(channel(mb.D, mb.B, mb.t, mb.T, mb.r, UKPFC_JH[mb.name]), { name: mb.name });
      return Object.assign(rhsHot(mb.D, mb.B, mb.t), { name: mb.name });
    }
    return { iSection, channel, rhsHot, angle, torsionOpen, classify, ownProps, UKPFC_JH };
  })();

  // ---------------------------------------------------------------- analysis engines
  function twistOpen(Tq, L, J, H, warping) {
    const { E, G } = materials, a = Math.sqrt(E * H / (G * J));
    let phi, p1, p2, p3;
    if (warping === 'Free') {
      const z1 = 0, z2 = L / 2, th = Math.tanh(L / (2 * a));
      p1 = Tq / (G * J) * a / L * (L * L / (2 * a) * (1 / L - 2 * z1 / (L * L)) + Math.sinh(z1 / a) - th * Math.cosh(z1 / a));
      p3 = Tq / (G * J * a * a) * a / L * (Math.sinh(z1 / a) - th * Math.cosh(z1 / a));
      phi = Tq * a / (G * J) * a / L * (L * L / (2 * a * a) * (z2 / L - z2 * z2 / (L * L)) + Math.cosh(z2 / a) - th * Math.sinh(z2 / a) - 1);
      p2 = Tq / (G * J * a) * a / L * (Math.cosh(z2 / a) - th * Math.sinh(z2 / a) - 1);
    } else {
      const z1 = 0, z2 = 0.2 * L, z3 = L / 2, k = (1 + Math.cosh(L / a)) / Math.sinh(L / a);
      p2 = Tq / (2 * G * J * a) * (k * Math.cosh(z1 / a) - 2 * a / L - Math.sinh(z1 / a));
      p3 = Tq / (2 * G * J * a * a) * (k * Math.sinh(z1 / a) - Math.cosh(z1 / a));
      p1 = Tq / (2 * G * J) * (k * Math.sinh(z2 / a) + (1 - 2 * z2 / L) - Math.cosh(z2 / a));
      phi = Tq * a / (2 * G * J) * (k * (Math.cosh(z3 / a) - 1) + z3 / a * (1 - z3 / L) - Math.sinh(z3 / a));
    }
    // signed values at the reference evaluation points (rad, rad/mm, rad/mm2, rad/mm3) are kept for the report
    return { a, phi: Math.abs(phi), phi1: Math.abs(p1), phi2: Math.abs(p2), phi3: Math.abs(p3), signed: { phi, p1, p2, p3 } };
  }
  const TABLE15 = [[1.0, Infinity], [1.25, 770], [1.5, 515], [1.75, 380], [2.0, 340], [2.5, 275], [3.0, 225], [4.0, 170]];
  function lambdaLimRHS(D, B, py) {
    const r = D / B;
    if (r <= 1) return Infinity;
    for (let i = 1; i < TABLE15.length; i++) {
      const [r0, l0] = TABLE15[i - 1], [r1, l1] = TABLE15[i];
      if (r <= r1) return (isFinite(l0) ? l0 + (l1 - l0) * (r - r0) / (r1 - r0) : l1) * 275 / py;
    }
    return TABLE15[TABLE15.length - 1][1] * 275 / py;
  }

  // ---------------------------------------------------------------- line loads (N/mm)
  // Support member: outer masonry leaf, or a slab edge (the reference "Slab": the unfactored slab reactions replace the
  // masonry and the additional loads are not used). Torsion beam: inner leaf, or a slab on the beam.
  function supportLoad(p, gG, gQ) {
    if (p.sb_load === 'Slab') return { uls: p.PG_ssb * gG + p.PQ_ssb * gQ, sls: p.PG_ssb + p.PQ_ssb, G: p.PG_ssb, Q: p.PQ_ssb };
    const G = p.b_msb * p.h_msb * p.rho_m_sb + p.PGadd_sb;
    return { uls: (p.b_msb * p.h_msb * p.rho_m_sb + p.PGadd_sb) * gG + p.PQadd_sb * gQ, sls: p.b_msb * p.h_msb * p.rho_m_sb + p.PGadd_sb + p.PQadd_sb, G, Q: p.PQadd_sb };
  }
  function beamLoad(p, gG, gQ) {
    if (p.mb_load === 'Slab') return { uls: p.PG_smb * gG + p.PQ_smb * gQ, sls: p.PG_smb + p.PQ_smb, G: p.PG_smb, Q: p.PQ_smb };
    const G = p.h_mmb * p.b_mmb * p.rho_m_mb + p.PGadd_mb;
    return { uls: (p.h_mmb * p.b_mmb * p.rho_m_mb + p.PGadd_mb) * gG + p.PQadd_mb * gQ, sls: p.h_mmb * p.b_mmb * p.rho_m_mb + p.PGadd_mb + p.PQadd_mb, G, Q: p.PQadd_mb };
  }
  // Direct loads (a page option, not a reference input): dead and imposed line loads typed instead of a leaf height and
  // density. They enter the model above as a leaf of zero height carrying the totals as the additional dead and
  // imposed load, so calc() itself is unchanged. which = 'sb' (support member) | 'mb' (torsion beam);
  // lists = {G: [N/mm, ...], Q: [N/mm, ...]} (1 kN/m = 1 N/mm). Returns a copy of p.
  function directApply(p, which, lists) {
    const sum = a => (a || []).reduce((t, v) => t + Number(v), 0);
    const q = Object.assign({}, p);
    q[which === 'sb' ? 'h_msb' : 'h_mmb'] = 0;
    q['PGadd_' + which] = sum(lists.G);
    q['PQadd_' + which] = sum(lists.Q);
    return q;
  }
  // Leaf height (mm) whose weight b × h × ρ equals the dead load G (N/mm): the equivalence the option is tested with.
  const directHeight = (G, b, rho) => G / (b * rho);
  // Angle soffit h_bf from the beam soffit: a negative value shortens the lever and the load arm (the reference, probed).
  const hbfOf = (p, torsion) => torsion && p.sb_method === 'Cantilever' ? Math.min(p.h_bf || 0, 0) : 0;

  // ---------------------------------------------------------------- orchestrator
  function calc(p) {
    const { E, G, RHO } = materials;
    const gG = p.gG, gQ = p.gQ, torsion = p.design_method === 'Masonry support and torsion';
    const c = p.c, emb = p.e_mb || 0, shim = p.shim || 0, o = {};
    const ld1 = supportLoad(p, gG, gQ), P1 = ld1.uls, P1SLS = ld1.sls;
    Object.assign(o, { P1, P1SLS });
    let lh, lv, tsb, r1, pysb, Asb, cyysb;
    if (p.sb_method === 'Plate') {
      lh = p.lh; tsb = p.tsb; pysb = materials.py(p.grade_sb, tsb); Asb = tsb * p.lplate; cyysb = lh / 2 - (p.lplate - lh) / 2;
    } else if (p.sb_method === 'Cantilever') {
      const long = (p.orientation || 'Long') === 'Long';
      lh = long ? p.angle.D : p.angle.B; lv = long ? p.angle.B : p.angle.D; tsb = p.angle.t;
      pysb = materials.py(p.grade_sb, tsb); Asb = p.angle.A; cyysb = p.angle.cyy;
    } else {
      lh = p.lh; lv = p.lv; tsb = p.tsb; r1 = p.r1sb; pysb = materials.stainless(p.grade_sb, p.gmss || 1, p.tsb);
      Asb = tsb * (lh + lv + PI_T * (tsb / 2 + r1) / 2);
      cyysb = (lv * tsb * tsb / 2 + lh * tsb * (lh - tsb) / 2) / ((lh + lv) * tsb - tsb * tsb);
    }
    Object.assign(o, { pysb, Asbu: Asb, cyysb });
    let mb, L, py, w2u, w2s, w3u, w3s, wULS, wSLS, Mx, Fv;
    if (torsion) {
      mb = p.mb; L = p.L;
      py = materials.py(p.grade_mb, mb.type === 'HR' ? mb.t : Math.max(mb.t, mb.T));
      const ld2 = beamLoad(p, gG, gQ);
      w2u = ld2.uls; w2s = ld2.sls;
      w3u = Asb * RHO * gG; w3s = Asb * RHO;
      wULS = P1 + w2u + w3u + mb.A * RHO * gG; wSLS = P1SLS + w2s + w3s + mb.A * RHO;
      Mx = wULS * L * L / 8 + p.Mx_add; Fv = wULS * L / 2 + p.Fv_add;
      Object.assign(o, { w2ULS: w2u, w3ULS: w3u, w4ULS: mb.A * RHO * gG, wULS, wSLS, Mx, Fv, py,
        w2SLS: w2s, w3SLS: w3s, w4SLS: mb.A * RHO });
    }
    // support member
    if (p.sb_method === 'Plate') {
      const e1 = p.e1, Zp = tsb * tsb / 6, Mc = 1.2 * Zp * pysb, Mxp = P1 * e1;
      // toe limit min((1 + dm/c) x 1 mm, 2 mm) as for angles (live reference run V1: dm 30, c 100 -> 1.30 mm)
      Object.assign(o, { Zp, dm: lh + emb - c, Mc_plate: Mc, Mx_plate: Mxp, alpha_ls: Mxp / Mc, delta_lim: c > 0 ? Math.min(1 + (lh + emb - c) / c, 2) : 2 });
      if (torsion) {
        const ye = (mb.D + tsb) * Asb / (2 * (mb.A + Asb));
        const Ixx = (mb.Ix + mb.A * ye * ye) + Asb * (mb.D / 2 + tsb / 2 - ye) ** 2, Zxx = Ixx / (mb.D / 2 + tsb - ye);
        const s1 = Mx / Zxx, rad = 4 * pysb * pysb - 3 * s1 * s1, cfp = rad > 0 ? Math.sqrt(rad) : pysb;
        Object.assign(o, { ye_all: ye, Ixx_all: Ixx, Zxx_all: Zxx, sigma1: s1, cfp, alpha_ts: Math.max((cfp * cfp - s1 * s1) / (2 * cfp * pysb), 0) });
      }
      const toe = (P, e) => { const al = e / lh; return al * al * (3 - al) / 6 * P * lh ** 3 / (E * tsb ** 3 / 12); };
      o.delta_toe = toe(P1SLS, e1);
      Object.assign(o, { am: e1, bm: lh - e1, al: e1 / lh, Ieff: tsb ** 3 / 12 });
      if (p.con_stage) {
        const e1c = c - emb + p.b_msb / 2 - shim, P1c = 1200 * p.b_msb * p.rho_m_sb * gG, P1cS = 1200 * p.b_msb * p.rho_m_sb;
        Object.assign(o, { e1c, P1c, P1cS, Mx_platec: P1c * e1c, alpha_lsc: P1c * e1c / Mc, delta_c: toe(P1cS, e1c),
          amc: e1c, bmc: lh - e1c, alc: e1c / lh });
      }
      const aw = p.s_weld / Math.SQRT2, pw = materials.pw(p.grade_sb, p.weld_class);
      const RA = P1 * Math.max(1 + 3 * e1 / (2 * mb.B / 2), 1.4), Rp = (lh + mb.B) * tsb * pysb, Rl = 2 * Rp / L;
      const Rh = P1 * e1 / (p.s_weld / 2 + tsb / 2), Rw = Math.hypot(RA, Rl, Rh);
      Object.assign(o, { aw, RA, Rp, Rl, Rh, Rweld: Rw, pweld: pw, pc_weld: aw * pw });
    } else {
      let l, am, bm, sigma, I, kB, ftMoment;
      o.dm = lh + shim + (torsion ? emb : 0) - c;
      if (p.sb_method === 'Cantilever') {
        const e1 = p.e1, lb = p.fixing === 'Weld top' ? lv : lv - p.eb, hb = hbfOf(p, torsion);
        l = lb + hb;
        const M1 = P1 * e1, Z = tsb * tsb / 6;
        am = e1 + hb; bm = lh - e1; ftMoment = M1;
        // toe limit min((1 + d_m / c) x 1 mm, 2 mm), as the reference prints for angles (worked example X03: 1.87 mm)
        Object.assign(o, { l, lb, hb, M1, Zp: Z, Mc_plate: 1.2 * pysb * Z, delta_lim: Math.min((1 + o.dm / c) * 1, 2) });
        if (p.fixing === 'Bolt') {
          const lam = p.Bbolt / (lh + lb), Bs = Math.min(2 * p.Bbolt / (0.4 + lam), p.Bbolt), Bd = Math.min(2 * p.Bbolt / (0.8 + lam), p.Bbolt);
          Object.assign(o, { lam_m: lam, Beff_str: Bs, Beff_def: Bd });
          sigma = p.Bbolt * M1 / (Bs * Z); I = Bd * tsb ** 3 / 12; kB = p.Bbolt;
        } else { sigma = M1 / Z; I = tsb ** 3 / 12; kB = 1; }
      } else {
        const lb = lv - p.eb; l = lb - (tsb + r1);
        const g = lv - lb, dm = o.dm;
        const e = Math.max(10, Math.min(p.b_msb / 2, Math.min(16.5 * tsb / (l + g) - 0.73, 1) * dm));
        am = c - emb - (shim + tsb + r1) + PI_T * (tsb / 2 + r1) / 2 + e; bm = dm - e;
        const lam = p.Bbolt / (l + am + bm), Bs = Math.min(2 * p.Bbolt / (0.4 + lam), p.Bbolt), Bd = Math.min(2 * p.Bbolt / (0.8 + lam), p.Bbolt);
        sigma = 6 * P1 * p.Bbolt * (c - emb + e - shim - 0.5 * tsb) / (Bs * tsb * tsb);
        I = Bd * tsb ** 3 / 12; kB = p.Bbolt; ftMoment = P1 * am;
        Object.assign(o, { lb, l, g, e, lam_m: lam, Beff_str: Bs, Beff_def: Bd, delta_lim: Math.min((1 + dm / c) * 1, 2) });
      }
      const pd = am / l, qd = bm / l;
      Object.assign(o, { am, bm, pdef: pd, qdef: qd, sigma_heel: sigma, Ieff: I,
        delta_toe: (kB * P1SLS * l ** 3 * pd) / (12 * E * I) * (4 * pd * pd + 6 * pd * qd + 3 * pd + 3 * qd) });
      if (p.fixing === 'Bolt') {
        const b = materials.bolt(p.bolt_grade, p.bolt_size), Fs = P1 * p.Bbolt, Ft = p.Bbolt * ftMoment / l;
        Object.assign(o, { As_bolt: b.As, Ps_bolt: b.Ps, Pt_bolt: b.Pt, Fs_bolt: Fs, Ft_bolt: Ft, F_comb: Fs / b.Ps + Ft / b.Pt });
      } else {
        const aw = p.s_weld / Math.SQRT2, f = Math.hypot(P1, ftMoment / l);
        Object.assign(o, { aw, f_weld: f, P_weld: f / (p.l_weld * aw), pweld: materials.pw(p.grade_sb, p.weld_class) });
      }
    }
    // support-member utilisations (limits as the reference applies them)
    if (p.sb_method === 'Plate') {
      o.util_sb = { heel: o.Mx_plate / o.Mc_plate, toe: o.delta_toe / o.delta_lim, weld: o.Rweld / o.pc_weld };
      if (torsion) o.util_sb.transverse = o.alpha_ts > 0 ? o.alpha_ls / o.alpha_ts : Infinity;
      if (p.con_stage) {
        Object.assign(o.util_sb, { heel_c: o.alpha_lsc, toe_c: o.delta_c / o.delta_lim });
        if (torsion) o.util_sb.transverse_c = o.alpha_ts > 0 ? o.alpha_lsc / o.alpha_ts : Infinity;
      }
    } else {
      o.util_sb = { heel: o.sigma_heel / o.pysb, toe: o.delta_toe / o.delta_lim };
      if (p.fixing === 'Bolt') Object.assign(o.util_sb, { bolt: Math.max(o.Fs_bolt / o.Ps_bolt, o.Ft_bolt / o.Pt_bolt), bolt_comb: o.F_comb / 1.4 });
      else o.util_sb.weld = o.P_weld / o.pweld;
    }
    if (!torsion) return o;
    // torque
    const tp = mb.type === 'HR' ? { e0: 0 } : sections.torsionOpen(mb);
    const t = mb.t;
    const ecc = mb.type === 'C'
      ? [c + (p.b_msb + t) / 2 - emb - tp.e0, -(p.b_mmb - t) / 2 - emb - tp.e0, cyysb - tp.e0 + t / 2]
      : [(mb.B + p.b_msb) / 2 + c - emb, (mb.B - p.b_mmb) / 2 - emb, mb.B / 2 + cyysb];
    // slab loads act at e1 from the heel (outside the beam face, plus any shim) and e2 inside the beam face
    const face = mb.type === 'C' ? t / 2 : mb.B / 2;
    if (p.sb_load === 'Slab') ecc[0] = p.e1 + shim + face - tp.e0;
    if (p.mb_load === 'Slab') ecc[1] = face - p.e2 - tp.e0;
    // Mirrored channel (toes towards the cavity; not a reference option): the shear centre moves inboard by
    // B - t + 2 e0, so every lever grows by that amount.
    if (mb.type === 'C' && p.mirror) { const sh = mb.B - t + 2 * tp.e0; for (let i = 0; i < 3; i++) ecc[i] += sh; o.mirrorShift = sh; }
    const TqULS = Math.abs(P1 * ecc[0] + w2u * ecc[1] + w3u * ecc[2]), TqSLS = Math.abs(P1SLS * ecc[0] + w2s * ecc[1] + w3s * ecc[2]);
    const Tq = TqULS * L, Tqu = TqSLS * L;
    Object.assign(o, { e0mb: tp.e0, e1mb: ecc[0], e2mb: ecc[1], e3mb: ecc[2], TqULS, TqSLS, Tq, Tqu });
    if (mb.type !== 'HR') Object.assign(o, { Wn0: tp.Wn0, Sw1: tp.Sw1, Qf: tp.Qf, Qw: tp.Qw });
    let tw;
    if (mb.type === 'HR') { const To = Tq / 2; tw = { To, Tav: To / 2, phi: To / 2 / (G * mb.J) * L / 2, phi1: 0, phi2: 0, phi3: 0 }; Object.assign(o, { To, Tav: To / 2 }); }
    else { tw = twistOpen(Tq, L, mb.J, mb.H, p.warping); o.a = tw.a; o.L_a = L / tw.a; o.twistSigned = tw.signed; }
    Object.assign(o, { phi_deg: tw.phi * DEG, phi1_deg: tw.phi1 * DEG * 1e3, phi2_deg: tw.phi2 * DEG * 1e6, phi3_deg: tw.phi3 * DEG * 1e9 });
    // BS 5950 member
    const cl = sections.classify(mb, py);
    if (cl.cls === 'Slender') {         // the reference: "Section is slender - beyond scope"; no member checks
      Object.assign(o, { Class: cl.cls, bt: cl.bt, dt: cl.dt, eps: Math.sqrt(275 / py), slender: true, util: { scope: Infinity } });
      return o;
    }
    const Av = mb.type === 'HR' ? mb.A * mb.D / (mb.D + mb.B) : mb.t * mb.D, Pv = 0.6 * py * Av;
    let Mcx = (cl.cls === 'Plastic' || cl.cls === 'Compact') ? Math.min(py * mb.Sx, 1.2 * py * mb.Zx) : py * mb.Zx;
    o.highShear = Fv > 0.6 * Pv;
    if (o.highShear) {
      // S_v of the shear area; RHS: A_v = A D / (D + B) taken as two webs, S_v = A_v^2 / (8 t) (as the reference)
      const rho = (2 * Fv / Pv - 1) ** 2, Sv = mb.type === 'HR' ? Av * Av / (8 * mb.t) : mb.t * mb.D * mb.D / 4;
      Object.assign(o, { rho_v: rho, Sv });
      Mcx = Math.min(py * (mb.Sx - rho * Sv), 1.2 * py * mb.Zx);
    }
    Object.assign(o, { Class: cl.cls, bt: cl.bt, dt: cl.dt, eps: Math.sqrt(275 / py), Av, Pv, Mcx });
    const LE = p.K_LT * p.LLT, lam = LE / mb.ry;
    let Mb;
    Object.assign(o, { LE, lam });
    if (mb.type === 'HR') o.lamLim = lambdaLimRHS(mb.D, mb.B, py);
    if (mb.type === 'HR' && lam <= o.lamLim) Mb = Mcx;
    else {
      // I and C: Annex B.2.3; RHS beyond the Table 15 limit: lambda_LT from Annex B.2.6.1 (as the reference, live run V6)
      const pl = cl.cls === 'Plastic' || cl.cls === 'Compact', bw = pl ? 1 : mb.Zx / mb.Sx;
      let v = null, lamLT;
      if (mb.type === 'HR') {
        const gb = (1 - mb.Iy / mb.Ix) * (1 - mb.J / (2.6 * mb.Ix)), phib = Math.sqrt(mb.Sx ** 2 * gb / (mb.A * mb.J));
        lamLT = 2.25 * Math.sqrt(phib * lam * bw);
        Object.assign(o, { gb, phib });
      } else {
        v = 1 / (1 + 0.05 * (lam / mb.x) ** 2) ** 0.25;
        lamLT = mb.u * v * lam * Math.sqrt(bw);
      }
      const lamL0 = 0.4 * Math.sqrt(PI_T ** 2 * E / py);
      const pE = PI_T ** 2 * E / lamLT ** 2, eta = Math.max(7 * (lamLT - lamL0) / 1000, 0), phiLT = (py + (eta + 1) * pE) / 2;
      const pb = pE * py / (phiLT + Math.sqrt(phiLT * phiLT - pE * py));
      Mb = pb * (pl ? mb.Sx : mb.Zx);
      Object.assign(o, { v, bw, lamLT, lamL0, pE, eta, phiLT, pb });
    }
    o.Mb = Mb;
    // SCI P057
    const amp = 1 + 0.5 * Mx * p.mLT / Mb, Myt = Mx * tw.phi, sbyt = Myt / mb.Zy, sbx = Mx / mb.Zx, pv = 0.6 * py;
    let ib, local, tmax;
    if (mb.type === 'HR') {
      ib = Mx * p.mLT / Mb + sbyt / py * amp; local = sbx + sbyt;
      const tbw = Fv * (mb.Sx / 2) / (mb.Ix * 2 * mb.t), tt = tw.To / mb.C, tvt = tt * amp;
      tmax = tbw + tvt;
      Object.assign(o, { tau_bw: tbw, tau_t: tt, tau_vt: tvt, tau: tmax });
    } else {
      const sw = E * tp.Wn0 * tw.phi2;
      ib = Mx * p.mLT / Mb + (sbyt + sw) / py * amp; local = sbx + sbyt + sw;
      const tbw = Fv * tp.Qw / (mb.Ix * t), tbf = Fv * tp.Qf / (mb.Ix * mb.T);
      const ttw = G * t * tw.phi1, ttf = G * mb.T * tw.phi1, twf = E * tp.Sw1 * tw.phi3 / mb.T;
      const tvw = ttw * amp, tvf = (ttf + twf) * amp;
      tmax = Math.max(tbw + tvw, tbf + tvf);
      Object.assign(o, { sig_w: sw, tau_bw: tbw, tau_bf: tbf, tau_tw: ttw, tau_tf: ttf, tau_wf: twf, tau_vtw: tvw, tau_vtf: tvf,
        tau_w: tbw + tvw, tau_f: tbf + tvf });
    }
    Object.assign(o, { amp, Myt, sig_byt: sbyt, sig_bx: sbx, ib, local, pv });
    // SLS
    const phis = tw.phi * Tqu / Tq;
    o.phi_sls_deg = phis * DEG;
    let twistUtil = o.phi_sls_deg / p.phi_lim;
    if (p.twist_check === 'Lever') { o.delta_h_sls = p.h_d * phis; twistUtil = o.delta_h_sls / p.d_h_lim; }
    o.dy = 5 * wSLS * L ** 4 / (384 * E * mb.Ix) + p.d_add;
    o.d_lim = Math.min(L / p.k_delta, p.d_lim_abs);
    o.util = { shear: Fv / Pv, moment: Mx / Mcx, LTB: Mx * p.mLT / Mb, buckling_ib: ib, local: local / py,
      shear_stress: tmax / pv, twist: twistUtil, deflection: o.dy / o.d_lim };
    if (p.twist_on === 0) delete o.util.twist;            // the reference "Twist check" switched off
    if (p.defl_on === 0) delete o.util.deflection;        // the reference "Deflection" check switched off
    return o;
  }

  // ---------------------------------------------------------------- Eurocode 3 (UK NA) support member
  // The support-member models are the ones above (the reference / SCI P110 / SCI P157); the resistances follow
  // EN 1993-1-1, EN 1993-1-4 and EN 1993-1-8 with the UK NA partial factors gamma_M0 = 1.0, gamma_M2 = 1.25.
  const ec3 = (() => {
    const GM0 = 1.0, GM2 = 1.25, E = 210000, E_SS = 200000;
    // EN 10025-2 Table 7: yield by thickness band; tensile strength for 3 <= t <= 100 mm (t < 3 mm: higher)
    const FY = { S275: [[16, 275], [40, 265], [63, 255], [80, 245], [100, 235], [150, 225]],
      S355: [[16, 355], [40, 345], [63, 335], [80, 325], [100, 315], [150, 295]] };
    const FU = { S275: [[3, 430], [100, 410], [150, 400]], S355: [[3, 510], [100, 470], [150, 450]] };
    const BETA_W = { S275: 0.85, S355: 0.9 };                 // EN 1993-1-8 Table 4.1
    // EN 1993-1-4 Table 2.1, austenitic: cold rolled strip (t <= 8) and hot rolled strip (8 < t <= 13.5)
    const SS = { '304': [[8, 230, 540], [13.5, 210, 520]], '304L': [[8, 220, 520], [13.5, 200, 520]],
      '316': [[8, 240, 530], [13.5, 220, 530]], '316L': [[8, 240, 530], [13.5, 220, 530]] };
    const band = (tab, t) => { for (const r of tab) if (t <= r[0]) return r; throw new Error(`t = ${t} mm is outside the tabulated range`); };
    const fy = (grade, t) => band(FY[grade], t)[1];
    const fu = (grade, t) => band(FU[grade], t)[1];
    const ssProps = (grade, t) => { const r = band(SS[grade], t); return { fy: r[1], fu: r[2] }; };
    const AS = { M12: 84.3, M16: 157, M20: 245, M24: 353 };
    const FUB = { '4.6': 400, '8.8': 800 };
    const AV = { '4.6': 0.6, '8.8': 0.6 };                     // EN 1993-1-8 Table 3.4, shear plane through the thread
    const d0 = d => d + (d <= 14 ? 1 : d <= 24 ? 2 : 3);       // EN 1090-2 normal clearance holes
    // Bolts in one row along the angle: shear, tension, combined and bearing on the angle leg (EN 1993-1-8 Table 3.4).
    function bolts(p, t, fyPly, fuPly, stainless, FvEd, FtEd) {
      const d = Number(p.bolt_size.slice(1)), As = AS[p.bolt_size], fub = FUB[p.bolt_grade], h = d0(d);
      const FvRd = AV[p.bolt_grade] * fub * As / GM2, FtRd = 0.9 * fub * As / GM2;
      const fuB = stainless ? Math.min(0.5 * fyPly + 0.6 * fuPly, fuPly) : fuPly;   // stainless: f_u,red [verify]
      const alphaB = Math.min(p.eb / (3 * h), fub / fuB, 1), k1 = Math.min(1.4 * p.Bbolt / h - 1.7, 2.5);
      const FbRd = k1 * alphaB * fuB * d * t / GM2;
      return { d, As, fub, d0: h, FvRd, FtRd, fuB, alphaB, k1, FbRd, FvEd, FtEd,
        u_shear: FvEd / FvRd, u_tension: FtEd / FtRd, u_comb: FvEd / FvRd + FtEd / (1.4 * FtRd), u_bearing: FvEd / FbRd };
    }
    // Fillet weld, simplified method (EN 1993-1-8 4.5.3.3): f_vw,d = f_u / (sqrt(3) beta_w gamma_M2), weaker part joined.
    function weld(s, parts) {
      const w = parts.reduce((m, q) => (!m || q.fu < m.fu ? q : m), null);
      const a = s / Math.SQRT2, fvwd = w.fu / (Math.sqrt(3) * w.betaW * GM2);
      return { a, fu: w.fu, betaW: w.betaW, grade: w.grade, fvwd, FwRd: a * fvwd };
    }
    // p: engine input with EN 1990 factors in gG / gQ; beam: {D, B, t, A, Ix, grade, Mx (N mm, ULS)}; opt: {gM0ss}
    function support(p, beam, opt) {
      opt = opt || {};
      const torsion = p.design_method === 'Masonry support and torsion';
      const c = p.c, emb = torsion ? (p.e_mb || 0) : 0, shim = p.shim || 0, o = { gM0: GM0, gM2: GM2 };
      const ld = supportLoad(p, p.gG, p.gQ), P1 = ld.uls, P1SLS = ld.sls;
      Object.assign(o, { P1, P1SLS });
      const toePlate = (P, e, lh, t, Em) => { const al = e / lh; return al * al * (3 - al) / 6 * P * lh ** 3 / (Em * t ** 3 / 12); };
      if (p.sb_method === 'Plate') {
        const t = p.tsb, lh = p.lh, e1 = p.e1, f_y = fy(p.grade_sb, t), fyd = f_y / GM0;
        const Wpl = t * t / 4, McRd = Wpl * fyd, MEd = P1 * e1;
        Object.assign(o, { fy: f_y, Wpl, McRd, MEd, dm: lh + emb - c, delta_lim: c > 0 ? Math.min(1 + (lh + emb - c) / c, 2) : 2, E });
        const Ap = t * p.lplate, ye = (beam.D + t) * Ap / (2 * (beam.A + Ap));
        const Ixx = beam.Ix + beam.A * ye * ye + Ap * (beam.D / 2 + t / 2 - ye) ** 2, Zxx = Ixx / (beam.D / 2 + t - ye);
        const s1 = beam.Mx / Zxx, rad = 4 * fyd * fyd - 3 * s1 * s1, cfp = rad > 0 ? Math.sqrt(rad) : fyd;
        const ats = Math.max((cfp * cfp - s1 * s1) / (2 * cfp * fyd), 0), als = MEd / McRd;
        Object.assign(o, { Ap, ye, Ixx, Zxx, sigma1: s1, cfp, alpha_ts: ats, alpha_ls: als,
          delta_toe: toePlate(P1SLS, e1, lh, t, E), am: e1, bm: lh - e1, al: e1 / lh, Ieff: t ** 3 / 12 });
        if (p.con_stage) {
          const e1c = c - emb + p.b_msb / 2 - shim, P1c = 1200 * p.b_msb * p.rho_m_sb * p.gG, P1cS = 1200 * p.b_msb * p.rho_m_sb;
          Object.assign(o, { e1c, P1c, P1cS, MEdc: P1c * e1c, alpha_lsc: P1c * e1c / McRd, delta_c: toePlate(P1cS, e1c, lh, t, E) });
        }
        const w = weld(p.s_weld, [{ grade: p.grade_sb, fu: fu(p.grade_sb, t), betaW: BETA_W[p.grade_sb] },
          { grade: beam.grade, fu: fu(beam.grade, beam.T), betaW: BETA_W[beam.grade] }]);
        const RA = P1 * Math.max(1 + 3 * e1 / (2 * beam.B / 2), 1.4), Rp = (lh + beam.B) * t * fyd, Rl = 2 * Rp / beam.L;
        const Rh = P1 * e1 / (p.s_weld / 2 + t / 2), Rw = Math.hypot(RA, Rl, Rh);
        Object.assign(o, { weld: w, RA, Rp, Rl, Rh, Rweld: Rw });
        o.util = { heel: MEd / McRd, toe: o.delta_toe / 2, weld: Rw / w.FwRd };
        if (torsion) o.util.transverse = ats > 0 ? als / ats : Infinity;
        if (p.con_stage) Object.assign(o.util, { heel_c: o.alpha_lsc, toe_c: o.delta_c / 2, transverse_c: ats > 0 ? o.alpha_lsc / ats : Infinity });
        return o;
      }
      // angles (P157 load model as the BS path; elastic stress at the heel)
      let t, lh, lv, l, am, bm, sigma, I, kB, ftMoment, fyd, Em, stainless = p.sb_method === 'SCI', fyPly, fuPly;
      if (!stainless) {
        const long = (p.orientation || 'Long') === 'Long';
        lh = long ? p.angle.D : p.angle.B; lv = long ? p.angle.B : p.angle.D; t = p.angle.t;
        fyPly = fy(p.grade_sb, t); fuPly = fu(p.grade_sb, t); fyd = fyPly / GM0; Em = E;
        const e1 = p.e1, M1 = P1 * e1, Z = t * t / 6, lb = p.fixing === 'Weld top' ? lv : lv - p.eb, hb = hbfOf(p, torsion);
        l = lb + hb; am = e1 + hb; bm = lh - e1; ftMoment = M1;
        Object.assign(o, { l, lb, hb, M1, Zel: Z, delta_lim: Math.min(1 + (lh + shim + emb - c) / c, 2), dm: lh + shim + emb - c });
        if (p.fixing === 'Bolt') {
          const lam = p.Bbolt / (lh + lb), Bs = Math.min(2 * p.Bbolt / (0.4 + lam), p.Bbolt), Bd = Math.min(2 * p.Bbolt / (0.8 + lam), p.Bbolt);
          Object.assign(o, { lam_m: lam, Beff_str: Bs, Beff_def: Bd });
          sigma = p.Bbolt * M1 / (Bs * Z); I = Bd * t ** 3 / 12; kB = p.Bbolt;
        } else { sigma = M1 / Z; I = t ** 3 / 12; kB = 1; }
      } else {
        lh = p.lh; lv = p.lv; t = p.tsb; const r1 = p.r1sb, mp = ssProps(p.grade_sb, t);
        fyPly = mp.fy; fuPly = mp.fu; o.gM0ss = opt.gM0ss || 1.1; fyd = fyPly / o.gM0ss; Em = E_SS;
        const lb = lv - p.eb; l = lb - (t + r1);
        const g = lv - lb, dm = lh + shim + emb - c;
        const e = Math.max(10, Math.min(p.b_msb / 2, Math.min(16.5 * t / (l + g) - 0.73, 1) * dm));
        am = c - emb - (shim + t + r1) + Math.PI * (t / 2 + r1) / 2 + e; bm = dm - e;
        const lam = p.Bbolt / (l + am + bm), Bs = Math.min(2 * p.Bbolt / (0.4 + lam), p.Bbolt), Bd = Math.min(2 * p.Bbolt / (0.8 + lam), p.Bbolt);
        sigma = 6 * P1 * p.Bbolt * (c - emb + e - shim - 0.5 * t) / (Bs * t * t);
        I = Bd * t ** 3 / 12; kB = p.Bbolt; ftMoment = P1 * am;
        Object.assign(o, { lb, l, g, e, dm, lam_m: lam, Beff_str: Bs, Beff_def: Bd, delta_lim: Math.min((1 + dm / c) * 1, 2) });
      }
      const pd = am / l, qd = bm / l;
      Object.assign(o, { fy: fyPly, fu: fuPly, fyd, E: Em, am, bm, pdef: pd, qdef: qd, sigma_heel: sigma, Ieff: I,
        delta_toe: (kB * P1SLS * l ** 3 * pd) / (12 * Em * I) * (4 * pd * pd + 6 * pd * qd + 3 * pd + 3 * qd) });
      o.util = { heel: sigma / fyd, toe: o.delta_toe / o.delta_lim };
      if (p.fixing === 'Bolt') {
        o.bolt = bolts(p, t, fyPly, fuPly, stainless, P1 * p.Bbolt, p.Bbolt * ftMoment / l);
        Object.assign(o.util, { bolt: Math.max(o.bolt.u_shear, o.bolt.u_tension), bolt_comb: o.bolt.u_comb, bearing: o.bolt.u_bearing });
      } else {
        const w = weld(p.s_weld, [{ grade: p.grade_sb, fu: fuPly, betaW: BETA_W[p.grade_sb] }].concat(
          torsion ? [{ grade: beam.grade, fu: fu(beam.grade, beam.T), betaW: BETA_W[beam.grade] }] : []));
        const f = Math.hypot(P1, ftMoment / l);
        Object.assign(o, { weld: w, f_weld: f, Fw_Ed: f / p.l_weld });
        o.util.weld = o.Fw_Ed / w.FwRd;
      }
      return o;
    }
    return { GM0, GM2, E, E_SS, fy, fu, BETA_W, ssProps, AS, FUB, d0, bolts, weld, support };
  })();

  const SMS = { materials, sections, calc, ec3, supportLoad, beamLoad, direct: { apply: directApply, height: directHeight }, version: '1.2.0' };
  if (typeof module !== 'undefined' && module.exports) module.exports = SMS;
  else root.SMS = SMS;
})(typeof window !== 'undefined' ? window : globalThis);
