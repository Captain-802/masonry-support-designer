/* Eurocode 3 (UK NA) path of the Masonry Support Designer.
   Torsion beam: the beam-v03 engine, loaded as classic scripts in the
   same global scope, driven through its state S: analyse() + checks() give EN 1993-1-1 member checks with
   SCI P385 torsion (bending + torsion cross-section, EN 1993-6 Annex A buckling interaction, FE eigenvalue M_cr
   with load heights). Support member: SMS.ec3.support (sms.js).
   Coordinates: d = distance from the beam's OUTER face (the cavity side), positive INWARD. beam-v03's own +x also
   points inward, so its eccentricity is e = d_load - d_shear-centre (mm) and z_g is + above the shear centre. */
(function (root) {
  'use strict';
  const FAM = { UKPFC: 'pfc', UKB: 'ub', UKC: 'uc', RHS: 'rhs', SHS: 'shs' };
  const PREFIX = { pfc: 'UKPFC ', ub: 'UKB ', uc: 'UKC ', rhs: 'RHS ', shs: 'SHS ' };
  const RHO_STEEL = 78.5e-6;                                   // N/mm3
  const norm = s => String(s).toLowerCase().replace(/\s+/g, '').replace(/\.0(?=x|$)/g, '');
  const keyOf = r => (Array.isArray(r) ? r[0] : r.key);
  function table(fam) {
    return { pfc: PFC, ub: UB, uc: UC, rhs: RHS, shs: SHS_HF }[FAM[fam]] || [];
  }
  // sections of one family as {key, name} in the library's own order
  function list(fam) {
    const f = FAM[fam];
    return table(fam).map(r => ({ key: keyOf(r), name: PREFIX[f] + keyOf(r).replace(/\s*x\s*/g, 'x') }));
  }
  function find(fam, name) {
    const l = list(fam), n = norm(String(name).replace(/^(UKPFC|UKB|UKC|RHS|SHS)\s+/i, ''));
    return l.find(x => norm(x.key) === n) || null;
  }
  function setSection(fam, key) {
    const f = FAM[fam];
    S.family = f;
    if (f === 'pfc') S.sectionKey = key;
    else if (f === 'ub') S.ubKey = key;
    else if (f === 'uc') S.ucKey = key;
    else if (f === 'shs') { S.shsKey = key; S.shsType = 'HF'; }
    else { S.rhsKey = key; S.rhsType = 'HF'; }
  }
  // beam-v03 section (cm-based library units) -> the designer's mm properties
  function beamProps(sec, name, grade) {
    const type = sec.kind === 'channel' ? 'C' : sec.kind === 'box' ? 'HR' : 'I';
    const e0 = type === 'C' ? (sec.tp && Number.isFinite(+sec.tp.e0) ? +sec.tp.e0 : 10 * sec.e0) : 0;
    return { type, name, grade, D: +sec.D, B: +sec.B, t: +sec.tw, T: +sec.tf, r: +sec.r || 0, A: sec.A * 100,
      Ix: sec.Ix * 1e4, Iy: sec.Iy * 1e4, e0, esc: sec.tp && sec.tp.esc != null ? +sec.tp.esc : null };
  }
  function shearCentre(mb, mirror) {
    if (mb.type !== 'C') return mb.B / 2;
    return mirror ? mb.B - mb.t / 2 + mb.e0 : mb.t / 2 - mb.e0;
  }
  /* p: the designer's engine input (N, mm) with EN 1990 factors in gG / gQ; st: {mb_family, key, zgInner:'top'|'sc'}.
     Returns beam-v03's analysis + checks and the EC3 support-member checks. */
  function run(p, st) {
    S = JSON.parse(JSON.stringify(DEMO));
    const mirror = !!(p.mirror && FAM[st.mb_family] === 'pfc');
    Object.assign(S, { code: 'EC3', grade: p.grade_mb, L: p.L / 1000, restraint: 'ltb', mcrMethod: 'eigen', memberName: 'Torsion beam',
      axial: 0, Mz: 0, hinges: [], LT: null, leFactor: null, destab: false, eccOn: true, pfcMirror: mirror, E: 210000,
      plate: { on: false, side: 'bottom', t: 10, outL: 0, outR: 0 }, divisor: p.k_delta, deflAbs: p.d_lim_abs });
    setSection(st.mb_family, st.key);
    const warp = p.warping === 'Fixed' ? { warp: 1 } : {};
    S.ends = endsPreset('ss', { e1: Object.assign({ ss: 100 }, warp), e2: Object.assign({ ss: 100 }, warp) });
    S.combos = [{ id: 'c1', label: `ULS: ${p.gG}G + ${p.gQ}Q (Eq 6.10)`, factors: { G: p.gG, Q: p.gQ, W: 0, E: 0 }, sls: false, on: true },
      { id: 's1', label: 'SLS: 1.0G + 1.0Q', factors: { G: 1, Q: 1, W: 0, E: 0 }, sls: true, on: true }];
    S.ltbRestraints = [];
    if (p.LLT < p.L - 1e-6) for (let x = p.LLT; x < p.L - 1e-6; x += p.LLT) S.ltbRestraints.push({ pos: +(x / 1000).toFixed(4), v: true, phi: true, vp: false, phip: false });
    const sec = activeSection(), mb = beamProps(sec, st.name, p.grade_mb), D = mb.D;
    const dsc = shearCentre(mb, mirror), emb = p.e_mb || 0, shim = p.fixing === 'Bolt' ? (p.shim || 0) : 0;
    // support member: bearing level of the outer leaf, self-weight and its centroid
    let zgOuter, Asb, dSb, zgSb;
    if (p.sb_method === 'Plate') { zgOuter = -D / 2; Asb = p.tsb * p.lplate; dSb = p.lplate / 2 - p.lh; zgSb = -D / 2 - p.tsb / 2; }
    else if (p.sb_method === 'Cantilever') {
      const a = p.angle, long = (p.orientation || 'Long') === 'Long';
      zgOuter = -D / 2 + a.t; Asb = a.A; dSb = -shim - (long ? a.cx : a.cyy); zgSb = -D / 2 + (long ? a.cyy : a.cx);
    } else {
      zgOuter = -D / 2 + p.tsb; Asb = p.tsb * (p.lh + p.lv + Math.PI * (p.tsb / 2 + p.r1sb) / 2);
      const cyy = (p.lv * p.tsb * p.tsb / 2 + p.lh * p.tsb * (p.lh - p.tsb) / 2) / ((p.lh + p.lv) * p.tsb - p.tsb * p.tsb);
      dSb = -shim - cyy; zgSb = -D / 2 + p.lv / 4;
    }
    // slab options: the slab edge on the support acts e1 (+ shim) outside the beam face, the slab on the
    // beam e2 inside it
    const slabOut = p.sb_load === 'Slab', slabIn = p.mb_load === 'Slab';
    const dOuter = slabOut ? -(p.e1 + shim) : -(p.c - emb + p.b_msb / 2), dInner = slabIn ? p.e2 : emb + p.b_mmb / 2;
    const zgInner = st.zgInner === 'sc' ? 0 : D / 2, Lm = p.L / 1000;
    const loads = [], add = (label, w, cs, d, zg) => { if (w > 0) loads.push({ type: 'udl', case: cs, label, x1: 0, x2: Lm, w, e: d - dsc, zg }); };
    const lo = SMS.supportLoad(p, 1, 1), li = SMS.beamLoad(p, 1, 1);
    add(slabOut ? 'slab edge on support' : 'outer leaf', lo.G, 'G', dOuter, zgOuter);
    add(slabOut ? 'slab edge on support, imposed' : 'outer leaf imposed', lo.Q, 'Q', dOuter, zgOuter);
    add(slabIn ? 'slab on beam' : 'inner leaf', li.G, 'G', dInner, zgInner);
    add(slabIn ? 'slab on beam, imposed' : 'inner leaf imposed', li.Q, 'Q', dInner, zgInner);
    add('support member', Asb * RHO_STEEL, 'G', dSb, zgSb);
    S.loads = loads;
    const a = analyse(), ch = checks(a);
    const swE = typeof selfWeightEccentricity === 'function' ? +selfWeightEccentricity(sec) || 0 : 0;
    const sup = SMS.ec3.support(p, { D, B: mb.B, t: mb.t, T: mb.T, A: mb.A, Ix: mb.Ix, grade: p.grade_mb, L: p.L, Mx: (ch.Mx || 0) * 1e6 },
      { gM0ss: st.gM0ss });
    return { S: JSON.parse(JSON.stringify(S)), sec, mb, a, ch, loads, dsc, swE, swKN: selfWeightValue(sec), sup, mirror };
  }
  // the same support checks for "support member only" (no beam)
  function supportOnly(p, st) { return SMS.ec3.support(p, null, { gM0ss: st.gM0ss }); }
  function brief(r) {
    try { return renderDesignBrief(r.a, r.ch, r.sec); } catch (e) { return '<p>beam-v03 report unavailable: ' + String(e.message || e) + '</p>'; }
  }
  root.EC3MS = { FAM, list, find, run, supportOnly, brief, beamProps, shearCentre };
})(typeof window !== 'undefined' ? window : globalThis);
