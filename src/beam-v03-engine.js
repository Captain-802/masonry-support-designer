/* beam-v03 engine (commit a1a525a) */

/* ==== beam-v03 module 1 ==== */
/* ===========================================================================
   1. COMPUTATION ENGINE  (validated against closed-form solutions)
   Units throughout the solver: mm, N, N mm, N/mm.  Up = positive.
   =========================================================================== */
function linsolve(A,b){
  // Both callers solve restrained, symmetric stiffness matrices. Diagonal
  // scaling removes the translation/rotation unit disparity; Cholesky then
  // detects mechanisms, including hinge layouts that pass a restraint count.
  const n=b.length, scale=A.map((r,i)=>Math.sqrt(r[i]));
  const fail=()=>{ throw new Error('Singular or unstable stiffness matrix (mechanism or numerically ill-conditioned layout). Check supports, hinges and closely spaced nodes.'); };
  if(scale.some(s=>!Number.isFinite(s)||s<=0)||b.some(v=>!Number.isFinite(v))) fail();
  const C=Array.from({length:n},()=>new Float64Array(n));
  for(let i=0;i<n;i++) for(let j=0;j<=i;j++){
    let v=A[i][j]/scale[i]/scale[j];
    for(let k=0;k<j;k++) v-=C[i][k]*C[j][k];
    if(i===j){ if(!Number.isFinite(v)||v<=1e-12) fail(); C[i][j]=Math.sqrt(v); }
    else C[i][j]=v/C[j][j];
  }
  const y=new Float64Array(n), z=new Float64Array(n);
  for(let i=0;i<n;i++){ let v=b[i]/scale[i]; for(let j=0;j<i;j++) v-=C[i][j]*y[j]; y[i]=v/C[i][i]; }
  for(let i=n-1;i>=0;i--){ let v=y[i]; for(let j=i+1;j<n;j++) v-=C[j][i]*z[j]; z[i]=v/C[i][i]; }
  return Array.from(z,(v,i)=>v/scale[i]);
}
function buildNodes(L,supports,loads,nSub,extra){
  const pts=new Set([0,L]);
  supports.forEach(s=>pts.add(+s.pos));
  (extra||[]).forEach(x=>{ if(Number.isFinite(+x)) pts.add(+x); });   // e.g. internal-hinge positions
  loads.forEach(ld=>{ if(ld.type==='point'||ld.type==='moment') pts.add(+ld.pos);
    else {pts.add(+ld.x1); pts.add(+ld.x2);} });
  const P=[...pts].filter(x=>x>=0&&x<=L).sort((a,b)=>a-b);
  const step=L/nSub; let nodes=[];
  for(let i=0;i<P.length-1;i++){ const a=P[i],b=P[i+1];
    const ns=Math.max(1,Math.ceil((b-a)/step));
    for(let k=0;k<ns;k++) nodes.push(a+(b-a)*k/ns); }
  nodes.push(P[P.length-1]);
  nodes=[...new Set(nodes.map(x=>+x.toFixed(6)))].sort((a,b)=>a-b);
  return nodes;
}
/* Restraint rank of the in-plane model (pure). A zero-energy displacement is
   piecewise linear across moment releases: w(x) = a + b x/L + sum c_h <x-h>/L.
   Every support contributes the rows it restrains: vertical translation
   (types 'pinned' and 'fixed') and rotation (types 'fixed' and 'guided' - a
   guided / sliding end holds the rotation and leaves the translation free).
   The model is stable when the rows span all (2 + number of releases)
   rigid-body freedoms. Returns {rank, needed, ok}. */
function beamRestraintRank(L,supports,hinges){
  const releases=[...new Set((hinges||[]).filter(x=>x>0&&x<L))].sort((a,b)=>a-b);
  const basis=x=>[1,x/L,...releases.map(h=>Math.max(0,(x-h)/L))];
  const rows=[];
  supports.forEach(s=>{
    if(s.type==='pinned'||s.type==='fixed') rows.push(basis(s.pos));
    if(s.type==='fixed'||s.type==='guided') rows.push([0,1,...releases.map(h=>s.pos>h?1:0)]);
  });
  let rank=0;
  for(let col=0;col<releases.length+2;col++){
    let pivot=rank;
    for(let j=rank;j<rows.length;j++) if(Math.abs(rows[j][col])>Math.abs(rows[pivot][col])) pivot=j;
    if(pivot>=rows.length||Math.abs(rows[pivot][col])<1e-10) continue;
    [rows[rank],rows[pivot]]=[rows[pivot],rows[rank]];
    const div=rows[rank][col]; rows[rank]=rows[rank].map(v=>v/div);
    for(let j=rank+1;j<rows.length;j++){ const v=rows[j][col]; for(let k=col;k<rows[j].length;k++) rows[j][k]-=v*rows[rank][k]; }
    rank++;
  }
  return {rank,needed:releases.length+2,ok:rank>=releases.length+2,releases};
}
function solveBeam(L,EI,supports,loads,nSub=120,hinges){
  // Support types: 'pinned' (vertical translation held), 'fixed' (translation +
  // rotation held), 'guided' (rotation held, translation free - a sliding /
  // guided end: reaction moment, no vertical reaction). A free end is simply
  // absent from the list. Check the restraint matrix directly: finite
  // round-off in a large FE matrix can otherwise disguise an exact mechanism
  // as a very flexible stable beam.
  supports.forEach(s=>{ if(!['pinned','fixed','guided'].includes(s.type)) throw new Error('solveBeam: unknown support type "'+s.type+'" (pinned | fixed | guided)'); });
  if(!beamRestraintRank(L,supports,hinges).ok) throw new Error('Unstable beam mechanism: the end restraints do not hold every segment separated by internal hinges.');
  // Internal hinges = major-axis (in-plane) MOMENT RELEASES. At a hinge node the
  // two adjacent elements get INDEPENDENT rotation DOFs (sharing the translation),
  // so the transmitted bending moment is zero and the slope is discontinuous; the
  // BMD/deflection redistribute. With no hinges the DOF numbering is identical to
  // the classic 2-DOF-per-node scheme (wDof[i]=2i, thDof[i]=2i+1), so results are
  // byte-identical to before. Lateral/torsional continuity is NOT affected here
  // (that lives in the Mcr eigensolver, which reads the redistributed BMD).
  const hp=(hinges||[]).map(x=>+(+x).toFixed(6));
  const nodes=buildNodes(L,supports,loads,nSub,hp);
  const nN=nodes.length;
  const idx=new Map(nodes.map((x,i)=>[+x.toFixed(6),i]));
  const hingeNode=new Set();
  hp.forEach(x=>{ const i=idx.get(+x.toFixed(6)); if(i!=null && i>0 && i<nN-1) hingeNode.add(i); });
  // DOF layout: one translation per node; one shared rotation per node, except
  // hinge nodes which get a left rotation and a right rotation.
  const wDof=new Array(nN), thDof=new Array(nN), thL=new Array(nN), thR=new Array(nN);
  let nd=0;
  for(let i=0;i<nN;i++){ wDof[i]=nd++; if(hingeNode.has(i)){ thL[i]=nd++; thR[i]=nd++; } else { thDof[i]=nd++; } }
  const ndof=nd;
  const rotR=i=> hingeNode.has(i)? thR[i] : thDof[i];   // rotation seen by the element to the RIGHT of node i
  const rotL=i=> hingeNode.has(i)? thL[i] : thDof[i];   // rotation seen by the element to the LEFT  of node i
  const K=Array.from({length:ndof},()=>new Array(ndof).fill(0));
  const F=new Array(ndof).fill(0);
  for(let e=0;e<nN-1;e++){
    const Le=nodes[e+1]-nodes[e], c=EI/(Le*Le*Le);
    const k=[[12,6*Le,-12,6*Le],[6*Le,4*Le*Le,-6*Le,2*Le*Le],
             [-12,-6*Le,12,-6*Le],[6*Le,2*Le*Le,-6*Le,4*Le*Le]];
    const d=[wDof[e], rotR(e), wDof[e+1], rotL(e+1)];
    for(let a=0;a<4;a++) for(let b=0;b<4;b++) K[d[a]][d[b]]+=k[a][b]*c;
  }
  loads.forEach(ld=>{
    if(ld.type==='point'){ const i=idx.get(+(+ld.pos).toFixed(6)); if(i!=null) F[wDof[i]]+=ld.P; }
    else if(ld.type==='moment'){ const i=idx.get(+(+ld.pos).toFixed(6)); if(i!=null) F[rotR(i)]+=ld.M; }
    else if(ld.type==='udl'){
      for(let e=0;e<nN-1;e++){ const xa=nodes[e],xb=nodes[e+1];
        if(xb<=ld.x1+1e-9||xa>=ld.x2-1e-9) continue; const Le=xb-xa;
        const wv=x=>{ if(ld.x2===ld.x1) return ld.w1; const t=(x-ld.x1)/(ld.x2-ld.x1); return ld.w1+(ld.w2-ld.w1)*t; };
        const wa=wv(xa),wb=wv(xb), d=[wDof[e], rotR(e), wDof[e+1], rotL(e+1)];
        F[d[0]]+=Le*(7*wa+3*wb)/20; F[d[1]]+=Le*Le*(3*wa+2*wb)/60;
        F[d[2]]+=Le*(3*wa+7*wb)/20; F[d[3]]+=-Le*Le*(2*wa+3*wb)/60;
      }
    }
  });
  // restrained DOFs per support type: the rotation DOF of a member end is the
  // node's single rotation (an end node is never a hinge node)
  const holdsW=t=>t==='pinned'||t==='fixed', holdsTh=t=>t==='fixed'||t==='guided';
  const fixed=new Set();
  supports.forEach(s=>{ const i=idx.get(+(+s.pos).toFixed(6)); if(i==null) return;
    if(holdsW(s.type)) fixed.add(wDof[i]); if(holdsTh(s.type)) fixed.add(rotR(i)); });
  const free=[]; for(let d=0;d<ndof;d++) if(!fixed.has(d)) free.push(d);
  const Kff=free.map(r=>free.map(c=>K[r][c])), Ff=free.map(r=>F[r]);
  const df=linsolve(Kff,Ff), d=new Array(ndof).fill(0);
  free.forEach((dof,j)=>d[dof]=df[j]);
  const R=new Array(ndof).fill(0);
  for(let i=0;i<ndof;i++){ let s=0; for(let j=0;j<ndof;j++) s+=K[i][j]*d[j]; R[i]=s-F[i]; }
  const reactions=supports.map(s=>{ const i=idx.get(+(+s.pos).toFixed(6));
    return {pos:s.pos,type:s.type,end:s.end,V:(holdsW(s.type)&&i!=null)?R[wDof[i]]:0,M:(holdsTh(s.type)&&i!=null)?R[rotR(i)]:0}; });
  const w=nodes.map((_,i)=>d[wDof[i]]);
  return {nodes,w,reactions};
}
/* Reaction moment of an end (fixed or guided) as the beam's end bending moment
   in the diagram convention (sagging positive, hogging negative), kN.m (pure).
   The nodal reaction R[rot] is anticlockwise-positive on the beam, so it is
   the negative of the beam moment at End 1 (x = 0) and equal to it at End 2
   (x = L): a fixed-fixed UDL prints -wL^2/12 at both ends, a guided End 2
   under a fixed - guided UDL prints +wL^2/6 (sagging). */
function reactionEndMomentKNm(r){ return ((r.end===2? r.M : -r.M)||0)/1e6; }
function sfdBmd(L,supports,loads,reactions,N=1000){
  const PF=[],PM=[];
  reactions.forEach(r=>{ PF.push([r.pos,r.V]); if(r.type==='fixed'||r.type==='guided') PM.push([r.pos,-r.M]); });
  // Applied point moments: with the solver's reaction convention (verified: R1=+M/L up,
  // V=+M/L), a ccw-positive applied moment REDUCES the sagging BM as x crosses it -
  // the diagram must close to zero at a pin. The previous +ld.M sign left the BMD
  // non-closing (M(L) = M instead of 0) for every applied-moment case; fixed-support
  // reaction moments (-r.M above) were and remain correct (fixed-fixed closures verified).
  loads.forEach(ld=>{ if(ld.type==='point') PF.push([ld.pos,ld.P]); else if(ld.type==='moment') PM.push([ld.pos,-ld.M]); });
  const crit=new Set([0,L]);
  PF.forEach(p=>{crit.add(p[0]-1e-4);crit.add(p[0]+1e-4);});
  PM.forEach(p=>{crit.add(p[0]-1e-4);crit.add(p[0]+1e-4);});
  loads.forEach(ld=>{ if(ld.type==='udl'){crit.add(ld.x1);crit.add(ld.x2);} });
  let grid=[]; for(let i=0;i<=N;i++) grid.push(L*i/N);
  grid=grid.concat([...crit].filter(x=>x>=0&&x<=L));
  grid=[...new Set(grid.map(x=>+x.toFixed(4)))].sort((a,b)=>a-b);
  const xs=[],V=[],M=[];
  grid.forEach(x=>{
    let v=0,m=0;
    PF.forEach(([p,f])=>{ if(p<=x+1e-9){ v+=f; m+=f*(x-p); } });
    PM.forEach(([p,mo])=>{ if(p<=x+1e-9) m+=mo; });
    loads.forEach(ld=>{ if(ld.type==='udl'){ const a=ld.x1,b=Math.min(ld.x2,x);
      if(b>a){ // exact closed-form integrals of the linear load (replaces the previous nn=40 trapezoid loop; makes triangular/trapezoidal BMD exact)
        const s=(ld.x2===ld.x1)?0:(ld.w2-ld.w1)/(ld.x2-ld.x1), u=b-a;
        const wA=(ld.x2===ld.x1)?ld.w1:ld.w1+s*(a-ld.x1);
        v+=wA*u+s*u*u/2;
        m+=wA*((x-a)*u-u*u/2)+s*((x-a)*u*u/2-u*u*u/3); } } });
    xs.push(x); V.push(v); M.push(m);
  });
  return {xs,V,M};
}
function pbFunc(lamLT,py,E){
  const lamL0=0.4*Math.sqrt(Math.PI*Math.PI*E/py);
  if(lamLT<=lamL0) return {pb:py,lamL0,pE:null,etaLT:0};
  const pE=Math.PI*Math.PI*E/(lamLT*lamLT);
  const etaLT=7*(lamLT-lamL0)/1000;
  const phiB=(py+(etaLT+1)*pE)/2;
  let pb=py*pE/(phiB+Math.sqrt(phiB*phiB-pE*py));
  return {pb:Math.min(pb,py),lamL0,pE,etaLT};
}
function pcFunc(lam,py,a,E){
  const lam0=0.2*Math.sqrt(Math.PI*Math.PI*E/py);
  if(lam<=lam0) return py;
  const pE=Math.PI*Math.PI*E/(lam*lam);
  const eta=a*(lam-lam0)/1000;
  const phi=(py+(eta+1)*pE)/2;
  return Math.min(py, py*pE/(phi+Math.sqrt(phi*phi-pE*py)));
}
function classify(bT,dt,eps,kind,opt){
  // BS 5950-1:2000 Table 11. Flange-outstand (rolled section) limit is common to
  // all open sections. The WEB limit differs by section type: a channel web has
  // its own flat 40e row for all three classes; an I/H web uses the "neutral
  // axis at mid-depth" row 80/100/120e for pure bending, or the "generally"
  // rows when axial COMPRESSION coexists (stress ratios per cl 3.5.5):
  //   Class 1: 80e/(1+r1),  Class 2: 100e/(1+1.5r1),  Class 3: 120e/(1+2r2),
  //   each but >= 40e;  r1 = Fc/(d t pyw) (-1 < r1 <= 1),  r2 = Fc/(Ag pyw).
  // At Fc = 0 these reduce exactly to 80/100/120e. Tension is conservatively
  // kept on the pure-bending row (relaxing limits for tension is not taken).
  const fc=bT<=9*eps?1:bT<=10*eps?2:bT<=15*eps?3:4;
  let wc, wlim=[80,100,120], webCase='bending', r1=null, r2=null;
  if(kind==='channel'){ wc = dt<=40*eps? 1 : 4; wlim=[40,40,40]; }
  else {
    if(opt && opt.Fc>0 && opt.dwt>0 && opt.Ag>0 && opt.py>0){
      r1=Math.min(opt.Fc/(opt.dwt*opt.py),1);
      r2=opt.Fc/(opt.Ag*opt.py);
      wlim=[Math.max(80/(1+r1),40), Math.max(100/(1+1.5*r1),40), Math.max(120/(1+2*r2),40)];
      webCase='bending+compression';
    }
    wc = dt<=wlim[0]*eps?1:dt<=wlim[1]*eps?2:dt<=wlim[2]*eps?3:4;
  }
  return {cls:Math.max(fc,wc),fc,wc,wlim,webCase,r1,r2};
}
function classifyBox(bt_flange,dt_web,eps,boxType){
  // BS 5950-1:2000 Table 12   box sections have SEPARATE flange (compression
  // due to bending, uses b) and web (neutral axis at mid-depth, uses d) rows.
  // For a square SHS the two ratios are numerically equal, so the flange row
  // (lower limits) always governs   this generalisation is fully backward
  // compatible with that case while correctly handling a true RHS (h != b).
  const flim = boxType==='CF'? [26,28,35]   : [28,32,40];
  const wlim = boxType==='CF'? [56,70,105]  : [64,80,120];
  const fc = bt_flange<=flim[0]*eps?1: bt_flange<=flim[1]*eps?2: bt_flange<=flim[2]*eps?3:4;
  const wc = dt_web<=wlim[0]*eps?1: dt_web<=wlim[1]*eps?2: dt_web<=wlim[2]*eps?3:4;
  return {cls:Math.max(fc,wc), fc, wc, flim, wlim};
}

/* ===========================================================================
   EUROCODE 3 (EN 1993-1-1 + UK National Annex) DESIGN FUNCTIONS
   Parallel to the BS5950 set above   same solver/diagrams feed both codes;
   only the cross-section/member checks differ. epsilon uses 235 (not 275).
   =========================================================================== */
function epsEC3(fy){ return Math.sqrt(235/fy); }

function classifyEC3(sec,eps,opt){
  // Table 5.2: box sections get flange="internal, compression" (33/38/42e) and
  // web="internal, bending" (72/83/124e). Open sections (I/H/channel) get
  // outstand flange (9/10/14e) and internal web, bending (72/83/124e)   EC3
  // does not carve out a separate channel-web row the way BS5950 Table 11 does.
  //
  // With coexistent axial COMPRESSION (opt.NEd in N, opt.fy in N/mm2) the web
  // must be classified as "internal part in bending and compression":
  //   Class 1: 396e/(13a-1) for a>0.5, 36e/a otherwise
  //   Class 2: 456e/(13a-1) for a>0.5, 41.5e/a otherwise
  //   Class 3: 42e/(0.67+0.33psi) for psi>-1
  // a from the plastic stress block with the web(s) conservatively carrying all
  // of N_Ed; psi from the conservative elastic form psi = 2N_Ed/(A f_y) - 1.
  // At N_Ed = 0 these reduce exactly to 72e / 82.9e / 124e, i.e. continuous
  // with the pure-bending row. Tension is conservatively ignored (kept at the
  // pure-bending limits, which are the tighter choice for a partly tensioned web).
  const flim = sec.isBox? [33,38,42] : [9,10,14];
  let wlim=[72,83,124], webCase='bending', alphaW=null, psiW=null;
  const NEd = (opt && opt.NEd>0)? opt.NEd : 0;
  if(NEd>0 && opt.fy>0){
    const twc = sec.tw||sec.tf;                 // box normalisation sets tw = t
    const nWeb = sec.isBox? 2 : 1;
    const cw = sec.dt*twc;                      // web flat depth c, mm
    alphaW = Math.min(0.5*(1 + NEd/(nWeb*cw*twc*opt.fy)), 1);
    psiW = Math.min(2*NEd/((sec.A*1e2)*opt.fy) - 1, 1);
    const l1 = alphaW>0.5? 396/(13*alphaW-1) : 36/alphaW;
    const l2 = alphaW>0.5? 456/(13*alphaW-1) : 41.5/alphaW;
    const l3 = psiW>-1? 42/(0.67+0.33*psiW) : 62*(1-psiW)*Math.sqrt(-psiW);
    wlim=[l1,l2,l3]; webCase='bending+compression';
  }
  let mzFlange=null;
  if(opt&&opt.minorBending){
    if(sec.kind==='I'){
      // I/H under M_z (19 Sep 2026 gap closure, item 1.9(b)): the web lies on the
      // z-z axis, so an applied M_z does not stress it - the web keeps its y-y
      // bending (+N) limits above. The flanges are outstands under the combined
      // stress: the outstand compressed by M_z has both its root and its tip in
      // compression (alpha = 1), for which Table 5.2 sheet 2 gives 9e/alpha,
      // 10e/alpha and 21e sqrt(k_sigma) with k_sigma >= 0.43 (EN 1993-1-5 Table
      // 4.2), i.e. never stricter than the uniform-compression bound 9e/10e/14e
      // already applied to flim; the opposite outstand (tip relieved by M_z) has
      // laxer limits and never governs. The uniform-compression bound is kept.
      mzFlange='outstand';
    } else {
      // Channels and hollow sections: a wall parallel to the web is compressed
      // over its full width by M_z - conservative uniform-compression bound.
      wlim=[33,38,42]; webCase='biaxial: uniform-compression web bound'; mzFlange='uniform-bound';
    }
  }
  const fc = sec.bT<=flim[0]*eps?1:sec.bT<=flim[1]*eps?2:sec.bT<=flim[2]*eps?3:4;
  const wc = sec.dt<=wlim[0]*eps?1:sec.dt<=wlim[1]*eps?2:sec.dt<=wlim[2]*eps?3:4;
  return {cls:Math.max(fc,wc),fc,wc,flim,wlim,webCase,alphaW,psiW,mzFlange};
}

function avEC3(sec){
  // EC3 cl 6.2.6(3) shear area, mm^2. Box: Av=A.D/(D+B), verified against the
  // verified commercial-software SHS example (collapses to A/2 for a square section, same as
  // the BS5950 rule). Channel/I-H formulas verified to ~98% against the PFC
  // example (residual gap likely a small area-rounding difference between the
  // BS5950-table and EC3-table figures for the same physical section).
  const A=sec.A*100;
  if(sec.isBox) return A*sec.D/(sec.D+sec.B);
  if(sec.kind==='channel') return A-2*sec.B*sec.tf+(sec.tw+sec.r)*sec.tf;
  const hw=sec.D-2*sec.tf;
  return Math.max(A-2*sec.B*sec.tf+(sec.tw+2*sec.r)*sec.tf, 1.2*hw*sec.tw);
}

// C1 for end-moment loading (linear gradient), psi = ratio of end moments.
// Source: SCI P360 Table 2.4 / section tables EC3 explanatory notes Table (8.1)  
// these are the precise published values (a fraction more refined than the
// older NCCI SN003 figures, e.g. ?=0.75?1.17 not 1.14; ?=-1.00?2.76, properly
// monotonic, not the 2.55 sometimes mis-cited).
const C1_END_MOMENT=[[1.00,1.00],[0.75,1.17],[0.50,1.36],[0.25,1.56],[0,1.77],[-0.25,2.00],[-0.50,2.24],[-0.75,2.49],[-1.00,2.76]];
// ---- NCCI SN006a-EN-EU: elastic critical moment of cantilevers ----
// Mcr = C * Mcr0;  Mcr0 = (pi/L)*sqrt(E*Iz*G*It);  kwt = (1/L)*sqrt(E*Iw/(G*It));
// eta = za/(hs/2), hs = h - tf (doubly symmetric I/H); +eta destabilising.
// C from Tables 3.1 (UDL q), 3.2 (tip point F), 3.3 (tip moment); warping at
// the support free or restrained. q+F combos by Eq (7). Values verbatim.
const SN006={
 eta:[-2,-1.5,-1,-0.5,0,0.25,0.5,0.75,1,1.25,1.5,2,2.5,3],
 kwt:[0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1],
 q:{free:[
  [2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04],
  [2.42,2.34,2.25,2.16,2.06,2.02,1.97,1.92,1.87,1.82,1.77,1.67,1.58,1.49],
  [2.87,2.71,2.53,2.34,2.13,2.03,1.92,1.82,1.71,1.61,1.52,1.35,1.20,1.07],
  [3.37,3.14,2.87,2.56,2.22,2.05,1.87,1.71,1.56,1.42,1.30,1.09,0.93,0.81],
  [3.93,3.62,3.25,2.82,2.32,2.06,1.82,1.60,1.41,1.25,1.12,0.91,0.76,0.65],
  [5.13,4.67,4.11,3.39,2.50,2.06,1.69,1.39,1.17,0.99,0.86,0.68,0.55,0.47],
  [6.40,5.79,5.04,4.02,2.66,2.02,1.54,1.21,0.98,0.81,0.70,0.54,0.43,0.36],
  [9.07,8.16,7.04,5.38,2.88,1.85,1.26,0.92,0.72,0.59,0.50,0.38,0.30,0.25],
  [11.8,10.6,9.11,6.82,3.00,1.65,1.03,0.73,0.56,0.46,0.38,0.29,0.23,0.19],
  [14.6,13.1,11.2,8.30,3.08,1.47,0.87,0.61,0.46,0.37,0.31,0.23,0.19,0.16]],
 restr:[
  [2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04],
  [2.81,2.71,2.61,2.50,2.39,2.33,2.27,2.21,2.15,2.09,2.04,1.92,1.81,1.70],
  [3.76,3.55,3.33,3.07,2.80,2.65,2.51,2.36,2.22,2.08,1.95,1.72,1.52,1.35],
  [4.80,4.50,4.15,3.72,3.23,2.97,2.71,2.47,2.23,2.03,1.84,1.54,1.31,1.13],
  [5.91,5.51,5.03,4.42,3.68,3.28,2.89,2.53,2.22,1.96,1.74,1.40,1.17,1.00],
  [8.22,7.63,6.91,5.93,4.57,3.82,3.14,2.59,2.16,1.84,1.59,1.24,1.02,0.86],
  [10.6,9.82,8.87,7.52,5.45,4.29,3.32,2.62,2.13,1.78,1.52,1.17,0.95,0.79],
  [15.5,14.3,12.9,10.8,7.15,5.10,3.65,2.74,2.17,1.78,1.51,1.15,0.92,0.77],
  [20.4,18.9,17.0,14.2,8.85,5.92,4.05,2.98,2.33,1.90,1.60,1.22,0.98,0.82],
  [25.4,23.5,21.1,17.6,10.6,6.80,4.54,3.30,2.57,2.09,1.76,1.33,1.07,0.89]]},
 F:{free:[
  [1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27],
  [1.39,1.37,1.34,1.31,1.28,1.26,1.25,1.23,1.21,1.19,1.17,1.13,1.08,1.04],
  [1.52,1.48,1.43,1.37,1.30,1.27,1.23,1.18,1.14,1.10,1.05,0.96,0.87,0.79],
  [1.65,1.60,1.53,1.45,1.34,1.27,1.21,1.13,1.06,0.99,0.92,0.80,0.69,0.61],
  [1.80,1.74,1.66,1.54,1.38,1.28,1.18,1.07,0.97,0.88,0.80,0.67,0.56,0.49],
  [2.15,2.07,1.94,1.75,1.45,1.27,1.10,0.94,0.81,0.70,0.61,0.49,0.40,0.34],
  [2.54,2.44,2.27,1.99,1.52,1.24,1.00,0.81,0.67,0.56,0.49,0.38,0.31,0.26],
  [3.41,3.26,3.01,2.52,1.60,1.14,0.81,0.61,0.49,0.40,0.34,0.26,0.21,0.17],
  [4.33,4.14,3.81,3.09,1.65,1.02,0.67,0.49,0.38,0.31,0.26,0.20,0.16,0.13],
  [5.29,5.06,4.63,3.70,1.68,0.91,0.57,0.40,0.31,0.25,0.21,0.16,0.13,0.11]],
 restr:[
  [1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.27],
  [1.55,1.52,1.49,1.45,1.42,1.40,1.37,1.35,1.33,1.31,1.28,1.24,1.19,1.14],
  [1.86,1.81,1.75,1.67,1.58,1.54,1.48,1.43,1.37,1.31,1.25,1.13,1.02,0.92],
  [2.20,2.13,2.04,1.93,1.77,1.68,1.58,1.47,1.36,1.26,1.16,0.99,0.85,0.74],
  [2.56,2.48,2.37,2.21,1.96,1.81,1.65,1.48,1.32,1.18,1.06,0.86,0.72,0.62],
  [3.36,3.26,3.10,2.82,2.35,2.03,1.72,1.44,1.21,1.04,0.90,0.71,0.58,0.49],
  [4.21,4.08,3.88,3.49,2.72,2.21,1.75,1.39,1.14,0.95,0.82,0.63,0.51,0.43],
  [5.99,5.82,5.52,4.90,3.46,2.53,1.84,1.39,1.10,0.91,0.77,0.59,0.47,0.40],
  [7.83,7.61,7.22,6.36,4.20,2.88,2.00,1.48,1.16,0.95,0.80,0.61,0.49,0.41],
  [9.69,9.43,8.94,7.84,4.98,3.27,2.21,1.62,1.26,1.03,0.86,0.66,0.53,0.44]]},
 M:{free:[0.50,0.50,0.50,0.51,0.51,0.52,0.53,0.54,0.54,0.54],
    restr:[0.50,0.50,0.50,0.51,0.52,0.55,0.59,0.68,0.80,0.93]}
};
function sn006C(tab,warp,kwt,eta){
  // bilinear interpolation; returns null outside the published ranges
  if(kwt<-1e-9||kwt>1+1e-9) return null;
  const K=SN006.kwt, T=SN006[tab][warp];
  let ki=0; while(ki<K.length-2 && kwt>K[ki+1]) ki++;
  const kf=(kwt-K[ki])/(K[ki+1]-K[ki]);
  if(tab==='M'){ return T[ki]+(T[ki+1]-T[ki])*Math.min(Math.max(kf,0),1); }
  const E2=SN006.eta;
  if(eta<E2[0]-1e-9||eta>E2[E2.length-1]+1e-9) return null;
  let ei=0; while(ei<E2.length-2 && eta>E2[ei+1]) ei++;
  const ef=(eta-E2[ei])/(E2[ei+1]-E2[ei]);
  const c1v=T[ki][ei]+(T[ki][ei+1]-T[ki][ei])*ef;
  const c2v=T[ki+1][ei]+(T[ki+1][ei+1]-T[ki+1][ei])*ef;
  return c1v+(c2v-c1v)*Math.min(Math.max(kf,0),1);
}
// ---- Serna et al. general C1 (SCI form, NSC Nov 2013, fork ends) ----
// C1 = sqrt( 35 Mmax^2 / (Mmax^2 + 9 M2^2 + 16 M3^2 + 9 M4^2) ).
// Note: the article's typeset equation omits the radical; restored here, as
// proven by its own worked examples (2.0 vs specialist software 2.03; 1.1 vs 1.1) and by
// the uniform-moment limit C1 = 1.
function sernaC1(Mmax,M2,M3,M4){
  const d=Mmax*Mmax+9*M2*M2+16*M3*M3+9*M4*M4;
  if(d<=0) return 1.0;
  return Math.sqrt(35*Mmax*Mmax/d);
}
function c1FromPsi(psi){
  psi=Math.max(-1,Math.min(1,psi));
  for(let i=0;i<C1_END_MOMENT.length-1;i++){
    const [p1,v1]=C1_END_MOMENT[i],[p2,v2]=C1_END_MOMENT[i+1];
    if(psi<=p1+1e-9 && psi>=p2-1e-9){ const t=(p1-psi)/(p1-p2); return v1+(v2-v1)*t; }
  }
  return 1.77;
}
function computeC1(a,isCant){
  // Returns {C1, method}. Uses the exact published SCI/section tables C1 values for
  // the pure end-moment-gradient case. For transverse loading or the genuinely
  // combined case   significant end moment AND transverse load together   SCI P360 itself states this
  // needs either the graphical NCCI SN003 method or software (specialist software), and
  // that "C1=1 is conservative" as the fallback; that is exactly what's used
  // here. Override below for a verified value (e.g. from specialist software, or a
  // closed-form fit such as Lopez/Yong/Serna that some commercial software uses).
  if(isCant) return {C1:1.0, method:'cantilever   safe default (the published tables assume a beam between two restrained ends)'};
  const Mmax=Math.abs(a.Mmax);
  if(Mmax<1e-9) return {C1:1.0, method:'negligible moment'};
  const endLevel=Math.max(Math.abs(a.M0end),Math.abs(a.MLend))/Mmax;
  if(endLevel<0.02){
    return {C1:1.0, method:'transverse load with no significant end moment   conservative default C1=1.0 used; override only with a verified value for the exact loading/restraint case'};
  }
  const psi=a.MLend/(a.M0end||1e-9);
  const C1endOnly=c1FromPsi(psi);
  return {C1:1.0, method:`combined end-moment (?=${psi.toFixed(2)}, end-moment-only C1 would be ${C1endOnly.toFixed(2)}) + transverse load   SCI P360/NCCI SN003 give this case only as graphs or via software (specialist software); "C1=1 is conservative" is SCI's own stated fallback, used here. Override below if you have a verified value (e.g. from specialist software, or a closed-form fit such as Lopez/Yong/Serna)`};
}

function mcrEC3(sec,Le,E,fy,C1){
  // Simplified (zg=0) Mcr for doubly symmetric open sections with load through
  // the shear centre. PFC/channel design is blocked elsewhere unless a
  // channel-specific Mcr method is implemented.
  const G=E/2.6;
  const Iz=sec.Iy*1e4, It=sec.J*1e4, Iw=(sec.Iw||0)*1e12;
  const term1=(Math.PI*Math.PI*E*Iz)/(Le*Le);
  const inside=Iw/Iz + (Le*Le*G*It)/(Math.PI*Math.PI*E*Iz);
  return C1*term1*Math.sqrt(Math.max(inside,0));
}
function ltbCurveEC3(sec){
  // Table 6.5 (cl 6.3.2.3, rolled/equivalent welded sections). Rolled I-section
  // (h/b<=2): curve b; (h/b>2): curve c. Anything not explicitly listed
  // (channels, box   though box is LTB-exempt anyway) defaults to curve d.
  if(sec.kind==='I'){ const hb=sec.D/sec.B; return hb<=2?{alphaLT:0.34,curve:'b'}:{alphaLT:0.49,curve:'c'}; }
  return {alphaLT:0.76,curve:'d'};
}
function chiLTEC3(lamLT,alphaLT,lamLT0,beta){
  const Phi=0.5*(1+alphaLT*(lamLT-lamLT0)+beta*lamLT*lamLT);
  return Math.min(1/(Phi+Math.sqrt(Math.max(Phi*Phi-beta*lamLT*lamLT,1e-12))),1.0);
}
function strutCurveEC3(sec,axis){
  // Table 6.2 (buckling curve allocation for axial compression).
  if(sec.isBox) return sec.boxType==='CF'? {alpha:0.49,curve:'c'} : {alpha:0.21,curve:'a'};
  if(sec.kind==='channel') return {alpha:0.49,curve:'c'};
  const hb=sec.D/sec.B;
  if(hb>1.2){ // rolled I-section
    if(sec.tf<=40) return axis==='y'?{alpha:0.21,curve:'a'}:{alpha:0.34,curve:'b'};
    return axis==='y'?{alpha:0.34,curve:'b'}:{alpha:0.49,curve:'c'};
  }
  // rolled H-section
  if(sec.tf<=100) return axis==='y'?{alpha:0.34,curve:'b'}:{alpha:0.49,curve:'c'};
  return {alpha:0.76,curve:'d'};
}
function chiStrutEC3(lamBar,alpha){
  const Phi=0.5*(1+alpha*(lamBar-0.2)+lamBar*lamBar);
  return Math.min(1/(Phi+Math.sqrt(Math.max(Phi*Phi-lamBar*lamBar,1e-12))),1.0);
}

function interpAt(xs,ys,xq){
  if(xq<=xs[0]) return ys[0];
  if(xq>=xs[xs.length-1]) return ys[ys.length-1];
  for(let i=0;i<xs.length-1;i++) if(xq>=xs[i]&&xq<=xs[i+1]){
    const t=(xs[i+1]===xs[i])?0:(xq-xs[i])/(xs[i+1]-xs[i]); return ys[i]+(ys[i+1]-ys[i])*t; }
  return ys[ys.length-1];
}
function mFactors(M2,M3,M4,MmaxSigned,M24abs){
  const Mm=Math.abs(MmaxSigned);
  if(Mm<1e-9) return {mLT:1,mx:1,a2:0,a3:0,a4:0};
  const s=MmaxSigned>=0?1:-1, a2=s*M2,a3=s*M3,a4=s*M4;
  const mLT=Math.max(0.2+(0.15*a2+0.5*a3+0.15*a4)/Mm,0.44);
  const mx =Math.max(0.2+(0.1*a2+0.6*a3+0.1*a4)/Mm,0.8*M24abs/Mm);
  return {mLT,mx,a2,a3,a4};
}


/* ==== beam-v03 module 2 ==== */
/* PFC channel section data and P385 torsion properties. */
const PFC=[
 // key,            mass,  D,   B,  tw,  tf,  r,  d,   bT,  dt,  Ix,   Iy,  rx,   ry,  Zx,  Zy,  Sx,  Sy,   u,    x,   e0,   Iw,       J,   A
 ["430x100x64",64.4,430,100,11,19,15,362,3.89,32.9,21900,722,16.3,2.97,1020,97.9,1220,176,0.916,22.5,3.27,0.219,63,82.1],
 ["380x100x54",54.0,380,100,9.5,17.5,15,315,4.31,33.2,15000,643,14.8,3.06,791,89.2,933,161,0.933,21.2,3.48,0.150,45.7,68.7],
 ["300x100x46",45.5,300,100,9,16.5,15,237,4.61,26.3,8230,568,11.9,3.13,549,81.7,641,148,0.944,17.0,3.67,0.0813,36.8,58.0],
 ["300x90x41",41.4,300,90,9,15.5,12,245,4.45,27.2,7220,404,11.7,2.77,481,63.1,568,114,0.934,18.4,3.18,0.0581,28.8,52.7],
 ["260x90x35",34.8,260,90,8,14,12,208,5.00,26.0,4730,353,10.3,2.82,364,56.3,425,102,0.943,17.2,3.31,0.0379,20.6,44.4],
 ["260x75x28",27.6,260,75,7,12,12,212,4.67,30.3,3620,185,10.1,2.30,278,34.4,328,62.0,0.932,20.5,2.61,0.0203,11.7,35.1],
 ["230x90x32",32.2,230,90,7.5,14,12,178,5.04,23.7,3520,334,9.27,2.86,306,55.0,355,98.9,0.949,15.1,3.45,0.0279,19.3,41.0],
 ["230x75x26",25.7,230,75,6.5,12.5,12,181,4.52,27.8,2750,181,9.17,2.35,239,34.8,278,63.2,0.945,17.3,2.77,0.0153,11.8,32.7],
 ["200x90x30",29.7,200,90,7,14,12,148,5.07,21.1,2520,314,8.16,2.88,252,53.4,291,94.5,0.952,12.9,3.60,0.0197,18.3,37.9],
 ["200x75x23",23.4,200,75,6,12.5,12,151,4.56,25.2,1960,170,8.11,2.39,196,33.8,227,60.6,0.956,14.7,2.91,0.0107,11.1,29.9],
 ["180x90x26",26.1,180,90,6.5,12.5,12,131,5.72,20.2,1820,277,7.40,2.89,202,47.4,232,83.5,0.950,12.8,3.63,0.0141,13.3,33.2],
 ["180x75x20",20.3,180,75,6,10.5,12,135,5.43,22.5,1370,146,7.27,2.38,152,28.8,176,51.8,0.945,15.3,2.85,0.00754,7.34,25.9],
 ["150x90x24",23.9,150,90,6.5,12,12,102,5.96,15.7,1160,253,6.18,2.89,155,44.4,179,76.9,0.937,10.8,3.71,0.00890,11.8,30.4],
 ["150x75x18",17.9,150,75,5.5,10,12,106,5.75,19.3,861,131,6.15,2.40,115,26.6,132,47.2,0.945,13.1,2.97,0.00467,6.1,22.8],
 ["125x65x15",14.8,125,65,5.5,9.5,12,82,5.00,14.9,483,80.0,5.07,2.06,77.3,18.8,89.9,33.2,0.942,11.1,2.54,0.00194,4.72,18.8],
 ["100x50x10",10.2,100,50,5,8.5,9,65,4.24,13.0,208,32.3,4.00,1.58,41.5,9.89,48.9,17.5,0.942,10.0,1.93,0.000491,2.53,13.0],
].map(a=>({key:a[0],mass:a[1],D:a[2],B:a[3],tw:a[4],tf:a[5],r:a[6],d:a[7],bT:a[8],dt:a[9],
  Ix:a[10],Iy:a[11],rx:a[12],ry:a[13],Zx:a[14],Zy:a[15],Sx:a[16],Sy:a[17],u:a[18],x:a[19],e0:a[20],Iw:a[21],J:a[22],A:a[23]}));
const PFCmap=Object.fromEntries(PFC.map(s=>[s.key,s]));

const TP385_PFC={"430x100x64":[65.7,0.916,0.219,127,67.2,780,561,491,32.7,53.4],"380x100x54":[47.6,0.893,0.15,109.6,63,608,407,347,34.8,57.9],"300x100x46":[38.4,0.729,0.0813,83.3,52.1,427,260,205,36.7,62.7],"300x90x41":[30,0.699,0.0581,76.5,45.2,336,219,178,31.8,53.3],"260x90x35":[21.4,0.669,0.0379,65,40.7,254,155,120,33.1,56.5],"260x75x28":[12.2,0.649,0.0203,56.3,32.3,161,108,90,26.1,43.6],"230x90x32":[20.1,0.59,0.0279,55.9,37.3,215,120,86,34.5,60],"230x75x26":[12.3,0.558,0.0153,47.9,30.1,140,84.6,67.3,27.7,47.4],"200x90x30":[19.1,0.508,0.0197,47,33.4,179,88.4,56.7,36,63.7],"200x75x23":[11.6,0.481,0.0107,40.3,27.2,115,62.5,46.3,29.1,50.9],"180x90x26":[13.8,0.506,0.0141,42.3,30.4,143,69.5,40.8,36.3,64.7],"180x75x20":[7.6,0.498,0.00754,36.8,24.2,89.3,50.8,35.8,28.5,49.6],"150x90x24":[12.2,0.426,0.0089,34.3,25.6,111.1,49.3,20.9,37.1,66.8],"150x75x18":[6.31,0.429,0.00467,29.8,20.8,67.9,34.8,21,29.7,52.8],"125x65x15":[4.89,0.313,0.00194,21.3,14.7,40.3,21.1,12.1,25.4,45.2],"100x50x10":[2.64,0.215,0.000491,12.9,8.83,16.9,8.99,5.21,19.3,34.1]};


/* ==== beam-v03 module 3 ==== */
/* SHS section data and P385 torsion properties. */
const SHS_HF=[
 // key, mass, D, t, A, I, r, Z, S, J, C, dt
 // Rows tagged "MS derived": key/D/t/mass from the reference software 2025 UK library (EN 10210 hot-finished, S355);
 // A, I, r, Z, S, J, C derived from the geometry with EN 10210-2:2006 Annex A (ro = 1.5t, ri = 1.0t; J and C
 // per the ctBoxEN10210() convention, mean corner radius 1.25t). Validated against every section tables row above:
 // max deviation 0.52 % (docs/SECTION_LIBRARY.md). dt = D/t - 3 (BS 5950 Table 11 flat width, as this file).
 ["20x20x2.0",1.1,20,2,1.397,0.7393,0.7275,0.7393,0.9299,1.218,1.067,7],  // MS derived
 ["20x20x2.5",1.32,20,2.5,1.683,0.8355,0.7046,0.8355,1.083,1.406,1.2,5],  // MS derived
 ["25x25x2.0",1.41,25,2,1.797,1.561,0.932,1.249,1.534,2.524,1.811,9.5],  // MS derived
 ["25x25x2.5",1.71,25,2.5,2.183,1.805,0.9093,1.444,1.816,2.973,2.084,7],  // MS derived
 ["25x25x3.0",2,25,3,2.543,1.999,0.8865,1.599,2.06,3.35,2.299,5.33],  // MS derived
 ["25x25x3.2",2.1,25,3.2,2.681,2.063,0.8773,1.65,2.147,3.481,2.37,4.81],  // MS derived
 ["30x30x2.5",2.11,30,2.5,2.683,3.329,1.114,2.219,2.737,5.4,3.215,9],  // MS derived
 ["30x30x3.0",2.47,30,3,3.143,3.743,1.091,2.495,3.138,6.164,3.601,7],  // MS derived
 ["30x30x3.2",2.61,30,3.2,3.321,3.888,1.082,2.592,3.286,6.442,3.736,6.38],  // MS derived
 ["40x40x2.5",2.89,40,2.5,3.683,8.538,1.523,4.269,5.141,13.59,6.223,13],  // MS derived
 ["40x40x2.9",3.35,40,2.9,4.213,9.539,1.505,4.77,5.809,15.32,6.931,10.8],  // MS derived
 ["40x40x3.0",3.41,40,3,4.34,9.78,1.5,4.89,5.97,15.7,7.1,10.3],
 ["40x40x3.2",3.61,40,3.2,4.6,10.2,1.49,5.11,6.28,16.5,7.42,9.5],
 ["40x40x3.6",4.01,40,3.6,5.103,11.07,1.473,5.536,6.881,18.07,8.008,8.11],  // MS derived; check availability
 ["40x40x4.0",4.39,40,4,5.59,11.8,1.45,5.91,7.44,19.5,8.54,7],
 ["40x40x5.0",5.28,40,5,6.73,13.4,1.41,6.68,8.66,22.5,9.6,5],
 ["40x40x5.6",5.93,40,5.6,7.369,14.06,1.381,7.032,9.284,23.97,10.08,4.14],  // MS derived; check availability
 ["40x40x6.3",6.52,40,6.3,8.067,14.68,1.349,7.339,9.905,25.36,10.5,3.35],  // MS derived
 ["50x50x2.5",3.68,50,2.5,4.683,17.46,1.931,6.986,8.295,27.46,10.23,17],  // MS derived
 ["50x50x2.9",4.26,50,2.9,5.373,19.67,1.913,7.868,9.427,31.17,11.49,14.2],  // MS derived
 ["50x50x3.0",4.35,50,3,5.54,20.2,1.91,8.08,9.7,32.1,11.8,13.7],
 ["50x50x3.2",4.62,50,3.2,5.88,21.2,1.9,8.49,10.2,33.8,12.4,12.6],
 ["50x50x3.6",5.14,50,3.6,6.543,23.17,1.882,9.269,11.28,37.2,13.47,10.9],  // MS derived; check availability
 ["50x50x4.0",5.64,50,4,7.19,25,1.86,9.99,12.3,40.4,14.5,9.5],
 ["50x50x5.0",6.85,50,5,8.73,28.9,1.82,11.6,14.5,47.6,16.7,7],
 ["50x50x5.6",7.69,50,5.6,9.609,30.83,1.791,12.33,15.74,51.32,17.76,5.93],  // MS derived; check availability
 ["50x50x6.0",7.99,50,6,10.17,31.98,1.773,12.79,16.48,53.6,18.39,5.33],  // MS derived
 ["50x50x6.3",8.31,50,6.3,10.6,32.8,1.76,13.1,17,55.2,18.8,4.94],
 ["50x50x7.1",9.37,50,7.1,11.64,34.53,1.722,13.81,18.28,58.95,19.8,4.04],  // MS derived; check availability
 ["50x50x8.0",10.3,50,8,12.75,36.01,1.68,14.4,19.5,62.32,20.6,3.25],  // MS derived; check availability
 ["60x60x2.9",5.17,60,2.9,6.533,35.22,2.322,11.74,13.91,55.3,17.2,17.7],  // MS derived
 ["60x60x3.0",5.29,60,3,6.74,36.2,2.32,12.1,14.3,56.9,17.7,17],
 ["60x60x3.2",5.62,60,3.2,7.16,38.2,2.31,12.7,15.2,60.2,18.6,15.8],
 ["60x60x3.6",6.27,60,3.6,7.983,41.88,2.291,13.96,16.76,66.49,20.37,13.7],  // MS derived; check availability
 ["60x60x4.0",6.9,60,4,8.79,45.4,2.27,15.1,18.3,72.5,22,12],
 ["60x60x5.0",8.42,60,5,10.7,53.3,2.23,17.8,21.9,86.4,25.7,9],
 ["60x60x5.6",9.45,60,5.6,11.85,57.38,2.201,19.13,23.87,93.93,27.64,7.71],  // MS derived; check availability
 ["60x60x6.0",9.87,60,6,12.57,59.89,2.182,19.96,25.11,98.63,28.81,7],  // MS derived
 ["60x60x6.3",10.3,60,6.3,13.1,61.6,2.17,20.5,26,102,29.6,6.52],
 ["60x60x7.1",11.6,60,7.1,14.48,65.84,2.132,21.95,28.22,110.2,31.56,5.45],  // MS derived; check availability
 ["60x60x8.0",12.5,60,8,16,69.7,2.09,23.2,30.4,118,33.4,4.5],
 ["70x70x3.0",6.24,70,3,7.943,59.02,2.726,16.86,19.87,92.19,24.76,20.3],  // MS derived; check availability
 ["70x70x3.2",6.63,70,3.2,8.441,62.31,2.717,17.8,21.04,97.58,26.11,18.9],  // MS derived
 ["70x70x3.6",7.4,70,3.6,9.42,68.6,2.7,19.6,23.3,108,28.7,16.4],
 ["70x70x4.0",8.15,70,4,10.39,74.69,2.681,21.34,25.54,118.2,31.16,14.5],  // MS derived
 ["70x70x5.0",9.99,70,5,12.7,88.5,2.64,25.3,30.8,142,36.8,11],
 ["70x70x5.6",11.2,70,5.6,14.09,95.94,2.61,27.41,33.68,155.2,39.75,9.5],  // MS derived; check availability
 ["70x70x6.0",11.8,70,6,14.97,100.6,2.591,28.73,35.53,163.5,41.6,8.67],  // MS derived
 ["70x70x6.3",12.3,70,6.3,15.6,104,2.58,29.7,36.9,169,42.9,8.11],
 ["70x70x7.1",13.8,70,7.1,17.32,111.9,2.542,31.97,40.28,184.5,46.13,6.86],  // MS derived; check availability
 ["70x70x8.0",15,70,8,19.2,120,2.5,34.2,43.8,200,49.2,5.75],
 ["70x70x8.8",16.6,70,8.8,20.71,125.7,2.464,35.92,46.6,211.8,51.6,4.95],  // MS derived; check availability
 ["80x80x3.0",7.18,80,3,9.143,89.82,3.134,22.46,26.3,139.6,33.04,23.7],  // MS derived
 ["80x80x3.2",7.63,80,3.2,9.721,94.95,3.125,23.74,27.88,147.9,34.9,22],  // MS derived
 ["80x80x3.6",8.53,80,3.6,10.9,105,3.11,26.2,31,164,38.5,19.2],
 ["80x80x4.0",9.41,80,4,12,114,3.09,28.6,34,180,41.9,17],
 ["80x80x5.0",11.6,80,5,14.7,137,3.05,34.2,41.1,217,49.8,13],
 ["80x80x5.6",13,80,5.6,16.33,148.8,3.018,37.19,45.17,238.4,54.09,11.3],  // MS derived; check availability
 ["80x80x6.0",13.6,80,6,17.37,156.4,3,39.1,47.76,251.8,56.78,10.3],  // MS derived
 ["80x80x6.3",14.2,80,6.3,18.1,162,2.99,40.5,49.7,262,58.7,9.7],
 ["80x80x7.1",16.1,80,7.1,20.16,175.6,2.951,43.89,54.47,286.2,63.5,8.27],  // MS derived; check availability
 ["80x80x8.0",17.5,80,8,22.4,189,2.91,47.3,59.5,312,68.3,7],
 ["80x80x8.8",19.4,80,8.8,24.23,200.1,2.873,50.02,63.66,332.4,72.04,6.09],  // MS derived; check availability
 ["80x80x10.0",21.6,80,10,26.93,213.9,2.818,53.47,69.3,360,76.82,5],  // MS derived
 ["90x90x3.6",9.66,90,3.6,12.3,152,3.52,33.8,39.7,237,49.7,22],
 ["90x90x4.0",10.7,90,4,13.6,166,3.5,37,43.6,260,54.2,19.5],
 ["90x90x5.0",13.1,90,5,16.7,200,3.45,44.4,53,316,64.8,15],
 ["90x90x5.6",14.7,90,5.6,18.57,218.1,3.427,48.46,58.34,346.9,70.65,13.1],  // MS derived; check availability
 ["90x90x6.0",15.5,90,6,19.77,229.8,3.409,51.07,61.79,367.1,74.34,12],  // MS derived
 ["90x90x6.3",16.2,90,6.3,20.7,238,3.4,53,64.3,382,77,11.3],
 ["90x90x7.1",18.3,90,7.1,23,259.6,3.36,57.7,70.79,419.5,83.69,9.68],  // MS derived; check availability
 ["90x90x8.0",20.1,90,8,25.6,281,3.32,62.6,77.6,459,90.5,8.25],
 ["90x90x8.8",22.1,90,8.8,27.75,299.1,3.283,66.46,83.36,491.5,95.95,7.23],  // MS derived; check availability
 ["90x90x10.0",24.7,90,10,30.93,322.3,3.228,71.61,91.27,536,103.1,6],  // MS derived
 ["100x100x3.6",10.8,100,3.6,13.74,211.7,3.924,42.33,49.49,328.5,62.33,24.8],  // MS derived; check availability
 ["100x100x4.0",11.9,100,4,15.2,232,3.91,46.4,54.4,361,68.2,22],
 ["100x100x5.0",14.7,100,5,18.7,279,3.86,55.9,66.4,439,81.8,17],
 ["100x100x5.6",16.5,100,5.6,20.81,306.1,3.836,61.23,73.19,484.1,89.45,14.9],  // MS derived; check availability
 ["100x100x6.0",17.4,100,6,22.17,323.2,3.818,64.64,77.61,513,94.3,13.7],  // MS derived
 ["100x100x6.3",18.2,100,6.3,23.2,336,3.8,67.1,80.9,534,97.8,12.9],
 ["100x100x7.1",20.5,100,7.1,25.84,367,3.768,73.4,89.24,588.6,106.7,11.1],  // MS derived; check availability
 ["100x100x8.0",22.6,100,8,28.8,400,3.73,79.9,98.2,646,116,9.5],
 ["100x100x8.8",24.9,100,8.8,31.27,426.2,3.692,85.24,105.7,694.4,123.4,8.36],  // MS derived; check availability
 ["100x100x10.0",27.4,100,10,34.9,462,3.64,92.4,116,761,133,7],
 ["120x120x4.0",14.4,120,4,18.39,410.3,4.723,68.38,79.71,635.1,100.8,27],  // MS derived; check availability
 ["120x120x5.0",17.8,120,5,22.7,498,4.68,83,97.6,777,122,21],
 ["120x120x5.6",20,120,5.6,25.29,547.4,4.653,91.24,107.9,858.1,133.8,18.4],  // MS derived; check availability
 ["120x120x6.0",21.2,120,6,26.97,579.4,4.635,96.57,114.7,911.2,141.4,17],  // MS derived
 ["120x120x6.3",22.2,120,6.3,28.2,603,4.62,100,120,950,147,16],
 ["120x120x7.1",25,120,7.1,31.52,662.9,4.586,110.5,132.5,1051,161.2,13.9],  // MS derived; check availability
 ["120x120x8.0",27.6,120,8,35.2,726,4.55,121,146,1160,176,12],
 ["120x120x8.8",30.4,120,8.8,38.31,779.1,4.51,129.8,158.3,1252,188.7,10.6],  // MS derived; check availability
 ["120x120x10.0",33.7,120,10,42.9,852,4.46,142,175,1380,206,9],
 ["120x120x11.0",37.2,120,11,46.66,907.6,4.41,151.3,188.4,1484,218.7,7.91],  // MS derived; check availability
 ["120x120x12.0",39.5,120,12,50.29,958.2,4.365,159.7,200.8,1578,230.5,7],  // MS derived
 ["120x120x12.5",40.9,120,12.5,52.1,982,4.34,164,207,1620,236,6.6],
 ["140x140x5.0",21,140,5,26.7,807,5.5,115,135,1250,170,25],
 ["140x140x5.6",23.5,140,5.6,29.77,890.5,5.469,127.2,149.4,1387,187,22],  // MS derived; check availability
 ["140x140x6.0",24.9,140,6,31.77,944.3,5.452,134.9,158.9,1475,198.1,20.3],  // MS derived
 ["140x140x6.3",26.1,140,6.3,33.3,984,5.44,141,166,1540,206,19.2],
 ["140x140x7.1",29.4,140,7.1,37.2,1086,5.403,155.1,184.4,1709,227.1,16.7],  // MS derived; check availability
 ["140x140x8.0",32.6,140,8,41.6,1200,5.36,171,204,1890,249,14.5],
 ["140x140x8.8",36,140,8.8,45.35,1287,5.327,183.8,221.5,2048,268,12.9],  // MS derived; check availability
 ["140x140x10.0",40,140,10,50.9,1420,5.27,202,246,2270,294,11],
 ["140x140x11.0",44.1,140,11,55.46,1516,5.228,216.6,265.6,2448,314.2,9.73],  // MS derived; check availability
 ["140x140x12.0",47,140,12,59.89,1609,5.183,229.9,284.3,2616,332.8,8.67],  // MS derived
 ["140x140x12.5",48.7,140,12.5,62.1,1650,5.16,236,293,2700,342,8.2],
 ["140x140x14.2",55.3,140,14.2,69.29,1790,5.083,255.8,322.2,2952,369,6.86],  // MS derived; check availability
 ["140x140x16.0",61.3,140,16,76.61,1916,5.001,273.7,350.2,3196,393.9,5.75],  // MS derived
 ["150x150x5.0",22.6,150,5,28.7,1000,5.9,134,156,1550,197,27],
 ["150x150x5.6",25.3,150,5.6,32.01,1106,5.878,147.5,172.6,1718,217,23.8],  // MS derived; check availability
 ["150x150x6.0",26.8,150,6,34.17,1174,5.86,156.5,183.7,1828,230,22],  // MS derived
 ["150x150x6.3",28.1,150,6.3,35.8,1220,5.85,163,192,1910,240,20.8],
 ["150x150x7.1",31.7,150,7.1,40.04,1352,5.811,180.3,213.5,2121,264.3,18.1],  // MS derived; check availability
 ["150x150x8.0",35.1,150,8,44.8,1490,5.77,199,237,2350,291,15.8],
 ["150x150x8.8",38.7,150,8.8,48.87,1608,5.736,214.4,257,2549,312.9,14],  // MS derived; check availability
 ["150x150x10.0",43.1,150,10,54.9,1770,5.68,236,286,2830,344,12],
 ["150x150x11.0",47.6,150,11,59.86,1902,5.637,253.6,309.2,3057,368.5,10.6],  // MS derived; check availability
 ["150x150x12.0",50.8,150,12,64.69,2023,5.592,269.7,331.4,3272,391.1,9.5],  // MS derived
 ["150x150x12.5",52.7,150,12.5,67.1,2080,5.57,277,342,3380,402,9],
 ["150x150x14.2",59.8,150,14.2,74.97,2262,5.492,301.5,376.9,3707,435.7,7.56],  // MS derived; check availability
 ["150x150x16.0",65.2,150,16,83.01,2430,5.41,324,410.7,4026,467,6.38],  // MS derived
 ["160x160x5.0",24.1,160,5,30.7,1230,6.31,153,178,1890,226,29],
 ["160x160x5.6",27,160,5.6,34.25,1353,6.286,169.2,197.6,2098,249.2,25.6],  // MS derived; check availability
 ["160x160x6.0",28.7,160,6,36.57,1437,6.268,179.6,210.4,2233,264.4,23.7],  // MS derived
 ["160x160x6.3",30.1,160,6.3,38.3,1500,6.26,187,220,2330,275,22.4],
 ["160x160x7.1",33.9,160,7.1,42.88,1659,6.22,207.4,244.7,2595,304.3,19.5],  // MS derived; check availability
 ["160x160x8.0",37.6,160,8,48,1830,6.18,229,272,2880,335,17],
 ["160x160x8.8",41.5,160,8.8,52.39,1978,6.144,247.2,295.2,3125,361.3,15.2],  // MS derived; check availability
 ["160x160x10.0",46.3,160,10,58.9,2190,6.09,273,329,3480,398,13],
 ["160x160x11.0",51,160,11,64.26,2349,6.046,293.6,356.1,3759,427.1,11.5],  // MS derived; check availability
 ["160x160x12.0",54.6,160,12,69.49,2502,6.001,312.8,382.1,4028,454.2,10.3],  // MS derived
 ["160x160x12.5",56.6,160,12.5,72.1,2580,5.98,322,395,4160,467,9.8],
 ["160x160x14.2",63.3,160,14.2,80.7,2810,5.9,351,436,4580,508,8.27],
 ["160x160x16.0",70.2,160,16,89.41,3028,5.82,378.5,476.1,4988,546.3,7],  // MS derived
 ["180x180x5.0",27.3,180,5,34.73,1765,7.129,196.1,227.3,2718,289.9,33],  // MS derived
 ["180x180x6.0",32.5,180,6,41.37,2077,7.085,230.8,269,3215,340.2,27],  // MS derived
 ["180x180x6.3",34,180,6.3,43.3,2170,7.07,241,281,3360,355,25.6],
 ["180x180x7.1",38.4,180,7.1,48.56,2404,7.037,267.2,313.5,3744,392.8,22.4],  // MS derived; check availability
 ["180x180x8.0",42.7,180,8,54.4,2660,7,296,349,4160,434,19.5],
 ["180x180x8.8",47,180,8.8,59.43,2880,6.961,320,379.5,4524,468.7,17.5],  // MS derived; check availability
 ["180x180x10.0",52.5,180,10,66.9,3190,6.91,355,424,5050,518,15],
 ["180x180x11.0",57.9,180,11,73.06,3441,6.863,382.4,459.7,5468,557.6,13.4],  // MS derived; check availability
 ["180x180x12.0",62.1,180,12,79.09,3677,6.818,408.5,494.3,5873,594.8,12],  // MS derived
 ["180x180x12.5",64.4,180,12.5,82.1,3790,6.8,421,511,6070,613,11.4],
 ["180x180x14.2",72.2,180,14.2,92,4150,6.72,462,566,6710,670,9.68],
 ["180x180x16.0",80.2,180,16,102,4500,6.64,500,621,7340,724,8.25],
 ["200x200x5.0",30.4,200,5,38.7,2450,7.95,245,283,3760,362,37],
 ["200x200x5.6",34.1,200,5.6,43.21,2710,7.92,271,314.1,4174,400.5,32.7],  // MS derived; check availability
 ["200x200x6.0",36.2,200,6,46.17,2883,7.902,288.3,334.9,4449,425.7,30.3],  // MS derived
 ["200x200x6.3",38,200,6.3,48.4,3010,7.89,301,350,4650,444,28.7],
 ["200x200x7.1",42.8,200,7.1,54.24,3345,7.853,334.5,390.9,5189,492.7,25.2],  // MS derived; check availability
 ["200x200x8.0",47.7,200,8,60.8,3710,7.81,371,436,5780,545,22],
 ["200x200x8.8",52.5,200,8.8,66.47,4021,7.778,402.1,474.3,6288,590.2,19.7],  // MS derived; check availability
 ["200x200x10.0",58.8,200,10,74.9,4470,7.72,447,531,7030,655,17],
 ["200x200x11.0",64.8,200,11,81.86,4829,7.68,482.9,576.5,7629,705.7,15.2],  // MS derived; check availability
 ["200x200x12.0",69.6,200,12,88.69,5171,7.635,517.1,620.9,8208,754.4,13.7],  // MS derived
 ["200x200x12.5",72.3,200,12.5,92.1,5340,7.61,534,643,8490,778,13],
 ["200x200x14.2",81.1,200,14.2,103,5870,7.54,587,714,9420,854,11.1],
 ["200x200x16.0",90.3,200,16,115,6390,7.46,639,785,10300,927,9.5],
 ["220x220x8.0",53,220,8,67.15,5002,8.63,454.7,531.8,7765,669.4,24.5],  // MS derived; check availability
 ["220x220x8.8",58.1,220,8.8,73.51,5430,8.595,493.7,579.7,8459,725.7,22],  // MS derived; check availability
 ["220x220x10.0",65.6,220,10,82.93,6050,8.542,550,649.8,9473,806.8,19],  // MS derived; check availability
 ["220x220x11.0",71.7,220,11,90.66,6546,8.497,595.1,706.6,10290,871.3,17],  // MS derived; check availability
 ["220x220x12.5",80.9,220,12.5,102.1,7254,8.43,659.5,789.1,11480,963.2,14.6],  // MS derived; check availability
 ["220x220x14.2",91,220,14.2,114.7,8007,8.354,727.9,878.6,12770,1060,12.5],  // MS derived; check availability
 ["220x220x16.0",102,220,16,127.8,8749,8.273,795.3,969,14050,1156,10.8],  // MS derived; check availability
 ["250x250x5.0",38.3,250,5,48.73,4861,9.987,388.9,446.8,7430,576.9,47],  // MS derived; check availability
 ["250x250x5.6",42.9,250,5.6,54.41,5399,9.961,431.9,497.6,8271,640,41.6],  // MS derived; check availability
 ["250x250x6.0",45.7,250,6,58.17,5752,9.943,460.1,531,8825,681.3,38.7],  // MS derived
 ["250x250x6.3",47.9,250,6.3,61,6010,9.93,481,556,9240,712,36.7],
 ["250x250x7.1",54,250,7.1,68.44,6701,9.895,536.1,621.6,10320,792.1,32.2],  // MS derived; check availability
 ["250x250x8.0",60.3,250,8,76.8,7460,9.86,596,694,11500,880,28.3],
 ["250x250x8.8",66.4,250,8.8,84.07,8107,9.82,648.6,757.6,12570,955.4,25.4],  // MS derived; check availability
 ["250x250x10.0",74.5,250,10,94.9,9060,9.77,724,851,14100,1070,22],
 ["250x250x11.0",82.1,250,11,103.9,9818,9.722,785.4,926.4,15350,1153,19.7],  // MS derived; check availability
 ["250x250x12.0",88.5,250,12,112.7,10560,9.678,844.4,1000,16570,1237,17.8],  // MS derived
 ["250x250x12.5",91.9,250,12.5,117,10900,9.66,873,1040,17200,1280,17],
 ["250x250x14.2",103,250,14.2,132,12100,9.58,967,1160,19100,1410,14.6],
 ["250x250x16.0",115,250,16,147,13300,9.5,1060,1280,21100,1550,12.6],
 ["260x260x6.3",49.9,260,6.3,63.5,6790,10.3,522,603,10400,773,38.3],
 ["260x260x7.1",56.2,260,7.1,71.28,7567,10.3,582.1,674.2,11650,860.5,33.6],  // MS derived; check availability
 ["260x260x8.0",62.8,260,8,80,8420,10.3,648,753,13000,956,29.5],
 ["260x260x8.8",69.1,260,8.8,87.59,9164,10.23,704.9,822.2,14190,1039,26.5],  // MS derived; check availability
 ["260x260x10.0",77.7,260,10,98.9,10200,10.2,788,924,15900,1160,23],
 ["260x260x11.0",85.5,260,11,108.3,11110,10.13,854.7,1006,17350,1255,20.6],  // MS derived; check availability
 ["260x260x12.5",95.8,260,12.5,122,12400,10.1,951,1130,19400,1390,17.8],
 ["260x260x14.2",108,260,14.2,137,13700,9.99,1060,1260,21700,1540,15.3],
 ["260x260x16.0",120,260,16,153,15100,9.91,1160,1390,23900,1690,13.3],
 ["300x300x6.0",55.1,300,6,70.17,10080,11.98,672,772.2,15410,996.9,47],  // MS derived
 ["300x300x6.3",57.8,300,6.3,73.6,10500,12,703,809,16100,1040,44.6],
 ["300x300x7.1",65.1,300,7.1,82.64,11780,11.94,785,905.6,18060,1163,39.3],  // MS derived; check availability
 ["300x300x8.0",72.8,300,8,92.8,13100,11.9,875,1010,20200,1290,34.5],
 ["300x300x8.8",80.2,300,8.8,101.7,14310,11.86,953.7,1107,22060,1409,31.1],  // MS derived; check availability
 ["300x300x10.0",90.2,300,10,115,16000,11.8,1070,1250,24800,1580,27],
 ["300x300x11.0",99.4,300,11,125.9,17420,11.76,1161,1359,27050,1710,24.3],  // MS derived; check availability
 ["300x300x12.0",107,300,12,136.7,18780,11.72,1252,1470,29250,1840,22],  // MS derived
 ["300x300x12.5",112,300,12.5,142,19400,11.7,1300,1530,30300,1900,21],
 ["300x300x14.2",126,300,14.2,160,21600,11.6,1440,1710,33900,2110,18.1],
 ["300x300x16.0",141,300,16,179,23900,11.5,1590,1900,37600,2330,15.8],
 ["350x350x8.0",85.4,350,8,109,21100,13.9,1210,1390,32400,1790,40.8],
 ["350x350x8.8",94,350,8.8,119.3,23060,13.9,1317,1522,35410,1950,36.8],  // MS derived; check availability
 ["350x350x10.0",106,350,10,135,25900,13.9,1480,1720,39900,2190,32],
 ["350x350x11.0",117,350,11,147.9,28180,13.81,1611,1874,43550,2376,28.8],  // MS derived; check availability
 ["350x350x12.0",126,350,12,160.7,30430,13.76,1739,2030,47150,2563,26.2],  // MS derived
 ["350x350x12.5",131,350,12.5,167,31500,13.7,1800,2110,48900,2650,25],
 ["350x350x14.2",148,350,14.2,189,35200,13.7,2010,2360,54900,2960,21.6],
 ["350x350x16.0",166,350,16,211,38900,13.6,2230,2630,61000,3260,18.9],
 ["350x350x19.0",190,350,19,247.7,44820,13.45,2561,3055,70760,3744,15.4],  // MS derived
 ["350x350x22.0",217,350,22,283.4,50270,13.32,2873,3460,80010,4187,12.9],  // MS derived
 ["350x350x25.0",242,350,25,318.3,55310,13.18,3161,3845,88750,4595,11],  // MS derived
 ["400x400x8.0",97.9,400,8,124.8,31860,15.98,1593,1830,48690,2363,47],  // MS derived
 ["400x400x8.8",108,400,8.8,136.9,34800,15.94,1740,2004,53290,2579,42.5],  // MS derived; check availability
 ["400x400x10.0",122,400,10,155,39100,15.9,1960,2260,60100,2900,37],
 ["400x400x11.0",134,400,11,169.9,42660,15.85,2133,2471,65670,3153,33.4],  // MS derived; check availability
 ["400x400x12.0",145,400,12,184.7,46130,15.8,2306,2679,71180,3405,30.3],  // MS derived
 ["400x400x12.5",151,400,12.5,192,47800,15.8,2390,2780,73900,3530,29],
 ["400x400x14.2",170,400,14.2,217,53500,15.7,2680,3130,83000,3940,25.2],
 ["400x400x16.0",191,400,16,243,59300,15.6,2970,3480,92400,4360,22],
 ["400x400x20.0",235,400,20,300,71500,15.4,3580,4250,112000,5240,17],
 ["400x400x22.0",251,400,22,327.4,77260,15.36,3863,4612,122100,5646,15.2],  // MS derived
 ["400x400x25.0",282,400,25,368.3,85380,15.23,4269,5141,135900,6223,13],  // MS derived
 ["450x450x12.0",162,450,12,208.7,66460,17.85,2954,3419,102200,4368,34.5],  // MS derived
 ["450x450x16.0",213,450,16,275,85860,17.67,3816,4459,133200,5620,25.1],  // MS derived
 ["450x450x19.0",250,450,19,323.7,99540,17.54,4424,5208,155400,6497,20.7],  // MS derived
 ["450x450x22.0",286,450,22,371.4,112500,17.4,5000,5929,176700,7324,17.5],  // MS derived
 ["450x450x25.0",321,450,25,418.3,124700,17.27,5544,6624,197200,8101,15],  // MS derived
 ["450x450x28.0",355,450,28,464.2,136300,17.14,6058,7292,216800,8832,13.1],  // MS derived
 ["450x450x32.0",399,450,32,524.1,150700,16.96,6696,8143,241700,9735,11.1],  // MS derived
 ["500x500x12.0",181,500,12,232.7,92030,19.89,3681,4248,141200,5451,38.7],  // MS derived
 ["500x500x16.0",238,500,16,307,119300,19.71,4771,5554,184400,7038,28.2],  // MS derived
 ["500x500x19.0",280,500,19,361.7,138600,19.58,5545,6498,215500,8159,23.3],  // MS derived
 ["500x500x22.0",320,500,22,415.4,157100,19.44,6283,7411,245600,9222,19.7],  // MS derived
 ["500x500x25.0",360,500,25,468.3,174600,19.31,6986,8295,274600,10230,17],  // MS derived
 ["500x500x28.0",399,500,28,520.2,191300,19.18,7653,9149,302600,11180,14.9],  // MS derived
 ["500x500x32.0",450,500,32,588.1,212300,19,8491,10240,338200,12370,12.6],  // MS derived
 ["500x500x36.0",498,500,36,654.3,231700,18.82,9269,11280,372000,13470,10.9],  // MS derived
 ["550x550x16.0",263,550,16,339,160400,21.75,5833,6769,247300,8616,31.4],  // MS derived
 ["550x550x19.0",309,550,19,399.7,186800,21.62,6793,7930,289500,10010,25.9],  // MS derived
 ["550x550x22.0",355,550,22,459.4,212100,21.49,7714,9058,330400,11340,22],  // MS derived
 ["550x550x25.0",399,550,25,518.3,236300,21.35,8594,10150,370100,12610,19],  // MS derived
 ["550x550x28.0",443,550,28,576.2,259500,21.22,9436,11220,408400,13810,16.6],  // MS derived
 ["550x550x32.0",500,550,32,652.1,288700,21.04,10500,12580,457500,15330,14.2],  // MS derived
 ["550x550x36.0",555,550,36,726.3,316100,20.86,11500,13890,504400,16740,12.3],  // MS derived
 ["550x550x40.0",608,550,40,798.8,341800,20.68,12430,15140,549000,18060,10.8],  // MS derived
 ["600x600x25.0",439,600,25,568.3,311100,23.4,10370,12200,485300,15230,21],  // MS derived
 ["600x600x28.0",487,600,28,632.2,342100,23.26,11400,13490,536300,16720,18.4],  // MS derived
 ["600x600x32.0",550,600,32,716.1,381600,23.08,12720,15160,601900,18600,15.8],  // MS derived
 ["600x600x36.0",611,600,36,798.3,418800,22.91,13960,16760,664900,20370,13.7],  // MS derived
 ["600x600x40.0",671,600,40,878.8,453900,22.73,15130,18310,725100,22030,12],  // MS derived
 ["700x700x25.0",517,700,25,668.3,504700,27.48,14420,16850,782900,21230,25],  // MS derived
 ["700x700x28.0",575,700,28,744.2,556600,27.35,15900,18670,867000,23380,22],  // MS derived
 ["700x700x32.0",651,700,32,844.1,623100,27.17,17800,21040,975800,26110,18.9],  // MS derived
 ["700x700x36.0",724,700,36,942.3,686500,26.99,19610,23330,1081000,28700,16.4],  // MS derived
 ["700x700x40.0",797,700,40,1039,746900,26.81,21340,25540,1182000,31160,14.5]  // MS derived
].map(a=>({key:a[0],mass:a[1],D:a[2],t:a[3],A:a[4],I:a[5],r:a[6],Z:a[7],S:a[8],J:a[9],C:a[10],dt:a[11]}));

const SHS_CF=[
 // key, mass, D, t, A, I, r, Z, S, J, C, dt
 // Rows tagged "MS derived": key/D/t/mass from the reference software 2025 UK library (EN 10219 cold-formed, S355);
 // A, I, r, Z, S, J, C derived with EN 10219-2:2006 Annex B (ro = 2.0t / 2.5t / 3.0t for t <= 6 / 6 < t <= 10 / t > 10,
 // ri = ro - t; J and C per the ctBoxEN10210() convention with mean corner radius (ro + ri)/2). Validated against
 // every section tables row above: max deviation 0.52 %. dt = D/t - 5 (BS 5950 cold-formed flat width, as this file).
 ["25x25x2.0",1.36,25,2,1.74,1.48,0.924,1.19,1.47,2.53,1.8,7.5],
 ["25x25x2.5",1.64,25,2.5,2.09,1.69,0.899,1.35,1.71,2.97,2.07,5],
 ["25x25x3.0",1.89,25,3,2.41,1.84,0.874,1.47,1.91,3.33,2.27,3.33],
 ["30x30x2.0",1.68,30,2,2.137,2.722,1.129,1.815,2.205,4.54,2.747,10],  // MS derived; check availability
 ["30x30x2.5",2.03,30,2.5,2.59,3.16,1.1,2.1,2.61,5.4,3.2,7],
 ["30x30x3.0",2.36,30,3,3.01,3.5,1.08,2.34,2.96,6.15,3.58,5],
 ["40x40x2.0",2.31,40,2,2.94,6.94,1.54,3.47,4.13,11.3,5.23,15],
 ["40x40x2.5",2.82,40,2.5,3.59,8.22,1.51,4.11,4.97,13.6,6.21,11],
 ["40x40x3.0",3.3,40,3,4.21,9.32,1.49,4.66,5.72,15.8,7.07,8.33],
 ["40x40x4.0",4.2,40,4,5.35,11.1,1.44,5.54,7.01,19.4,8.48,5],
 ["50x50x2.0",2.93,50,2,3.737,14.15,1.946,5.659,6.662,22.63,8.513,20],  // MS derived
 ["50x50x2.5",3.6,50,2.5,4.59,16.9,1.92,6.78,8.07,27.5,10.2,15],
 ["50x50x3.0",4.25,50,3,5.41,19.5,1.9,7.79,9.39,32.1,11.8,11.7],
 ["50x50x4.0",5.45,50,4,6.95,23.7,1.85,9.49,11.7,40.4,14.4,7.5],
 ["50x50x5.0",6.56,50,5,8.36,27,1.8,10.8,13.7,47.5,16.6,5],
 ["60x60x3.0",5.19,60,3,6.61,35.1,2.31,11.7,14,57.1,17.7,15],
 ["60x60x4.0",6.71,60,4,8.55,43.6,2.26,14.5,17.6,72.6,22,10],
 ["60x60x5.0",8.13,60,5,10.4,50.5,2.21,16.8,20.9,86.4,25.6,7],
 ["60x60x6.0",9.2,60,6,12.03,56.07,2.159,18.69,23.68,98.41,28.62,5],  // MS derived; check availability
 ["70x70x2.5",5.17,70,2.5,6.59,49.4,2.74,14.1,16.5,78.5,21.2,23],
 ["70x70x3.0",6.13,70,3,7.81,57.5,2.71,16.4,19.4,92.4,24.7,18.3],
 ["70x70x3.5",7.06,70,3.5,8.99,65.1,2.69,18.6,22.2,106,28,15],
 ["70x70x4.0",7.97,70,4,10.1,72.1,2.67,20.6,24.8,119,31.1,12.5],
 ["70x70x5.0",9.7,70,5,12.4,84.6,2.62,24.2,29.6,142,36.7,9],
 ["70x70x6.0",11.1,70,6,14.43,95.17,2.568,27.19,33.83,163.5,41.41,6.67],  // MS derived; check availability
 ["80x80x3.0",7.07,80,3,9.01,87.8,3.12,22,25.8,140,33,21.7],
 ["80x80x3.5",8.16,80,3.5,10.4,99.8,3.1,25,29.5,161,37.6,17.9],
 ["80x80x4.0",9.22,80,4,11.7,111,3.07,27.8,33.1,180,41.8,15],
 ["80x80x5.0",11.3,80,5,14.4,131,3.03,32.9,39.7,218,49.7,11],
 ["80x80x6.0",13.2,80,6,16.8,149,2.98,37.3,45.8,252,56.6,8.33],
 ["90x90x3.0",8.01,90,3,10.2,127,3.53,28.3,33,201,42.5,25],
 ["90x90x3.5",9.26,90,3.5,11.8,145,3.51,32.2,37.9,232,48.5,20.7],
 ["90x90x4.0",10.5,90,4,13.3,162,3.48,36,42.6,261,54.2,17.5],
 ["90x90x5.0",12.8,90,5,16.4,193,3.43,42.9,51.4,316,64.7,13],
 ["90x90x6.0",15.1,90,6,19.2,220,3.39,49,59.5,368,74.2,10],
 ["100x100x3.0",8.96,100,3,11.4,177,3.94,35.4,41.2,279,53.2,28.3],
 ["100x100x4.0",11.7,100,4,14.9,226,3.89,45.3,53.3,362,68.1,20],
 ["100x100x5.0",14.4,100,5,18.4,271,3.84,54.2,64.6,441,81.7,15],
 ["100x100x6.0",17,100,6,21.6,311,3.79,62.3,75.1,514,94.1,11.7],
 ["100x100x8.0",21.4,100,8,27.2,366,3.67,73.2,91.1,645,114,7.5],
 ["120x120x3.0",10.8,120,3,13.81,312.3,4.756,52.06,60.24,487.7,78.15,35],  // MS derived; check availability
 ["120x120x4.0",14.2,120,4,18.1,402,4.71,67,78.3,637,101,25],
 ["120x120x5.0",17.5,120,5,22.4,485,4.66,80.9,95.4,778,122,19],
 ["120x120x6.0",20.7,120,6,26.4,562,4.61,93.7,112,913,141,15],
 ["120x120x8.0",26.4,120,8,33.6,677,4.49,113,138,1160,175,10],
 ["120x120x10.0",31.8,120,10,40.6,777,4.38,129,162,1380,203,7],
 ["140x140x4.0",16.8,140,4,21.3,652,5.52,93.1,108,1020,140,30],
 ["140x140x5.0",20.7,140,5,26.4,791,5.48,113,132,1260,170,23],
 ["140x140x6.0",24.5,140,6,31.2,920,5.43,131,155,1480,198,18.3],
 ["140x140x8.0",31.4,140,8,40,1130,5.3,161,194,1900,248,12.5],
 ["140x140x10.0",38.1,140,10,48.6,1310,5.2,187,230,2270,291,9],
 ["150x150x4.0",18,150,4,22.9,808,5.93,108,125,1270,162,32.5],
 ["150x150x5.0",22.3,150,5,28.4,982,5.89,131,153,1550,197,25],
 ["150x150x6.0",26.4,150,6,33.6,1150,5.84,153,180,1830,230,20],
 ["150x150x8.0",33.9,150,8,43.2,1410,5.71,188,226,2360,289,13.8],
 ["150x150x10.0",41.3,150,10,52.6,1650,5.61,220,269,2840,341,10],
 ["160x160x4.0",19.3,160,4,24.5,987,6.34,123,143,1540,185,35],
 ["160x160x5.0",23.8,160,5,30.4,1200,6.29,150,175,1900,226,27],
 ["160x160x6.0",28.3,160,6,36,1410,6.25,176,206,2240,264,21.7],
 ["160x160x8.0",36.5,160,8,46.4,1740,6.12,218,260,2900,334,15],
 ["160x160x10.0",44.4,160,10,56.6,2050,6.02,256,311,3490,395,11],
 ["180x180x5.0",27,180,5,34.4,1740,7.11,193,224,2720,290,31],
 ["180x180x6.0",32.1,180,6,40.8,2040,7.06,226,264,3220,340,25],
 ["180x180x6.3",33.3,180,6.3,42.41,2096,7.03,232.8,273.1,3383,354.1,23.6],  // MS derived; check availability
 ["180x180x8.0",41.5,180,8,52.8,2550,6.94,283,336,4190,432,17.5],
 ["180x180x10.0",50.7,180,10,64.6,3020,6.84,335,404,5070,515,13],
 ["180x180x12.0",58.5,180,12,74.5,3320,6.68,369,454,5870,584,10],
 ["180x180x12.5",60.5,180,12.5,77,3410,6.65,378,467,6050,600,9.4],
 ["200x200x5.0",30.1,200,5,38.4,2410,7.93,241,279,3760,362,35],
 ["200x200x6.0",35.8,200,6,45.6,2830,7.88,283,330,4460,426,28.3],
 ["200x200x6.3",37.2,200,6.3,47.45,2922,7.847,292.2,341.2,4682,443.5,26.7],  // MS derived; check availability
 ["200x200x8.0",46.5,200,8,59.2,3570,7.76,357,421,5820,544,20],
 ["200x200x10.0",57,200,10,72.6,4250,7.65,425,508,7070,651,15],
 ["200x200x12.0",66,200,12,84.1,4730,7.5,473,576,8230,743,11.7],
 ["200x200x12.5",68.3,200,12.5,87,4860,7.47,486,594,8500,765,11],
 ["250x250x6.0",45.2,250,6,57.6,5670,9.92,454,524,8840,681,36.7],
 ["250x250x6.3",47.1,250,6.3,60.05,5873,9.889,469.8,544.4,9290,711.2,34.7],  // MS derived; check availability
 ["250x250x8.0",59.1,250,8,75.2,7230,9.8,578,676,11600,878,26.3],
 ["250x250x10.0",72.7,250,10,92.6,8710,9.7,697,822,14200,1060,20],
 ["250x250x12.0",84.8,250,12,108,9860,9.55,789,944,16700,1230,15.8],
 ["250x250x12.5",88,250,12.5,112,10200,9.52,813,975,17300,1270,15],
 ["300x300x6.0",54.7,300,6,69.6,9960,12,664,764,15400,997,45],
 ["300x300x6.3",57,300,6.3,72.65,10340,11.93,689.5,794.9,16220,1042,42.6],  // MS derived; check availability
 ["300x300x8.0",71.6,300,8,91.2,12800,11.8,853,991,20300,1290,32.5],
 ["300x300x10.0",88.4,300,10,113,15500,11.7,1040,1210,25000,1570,25],
 ["300x300x12.0",104,300,12,132,17800,11.6,1180,1400,29500,1830,20],
 ["300x300x12.5",108,300,12.5,137,18300,11.6,1220,1450,30600,1890,19],
 ["350x350x6.0",63.8,350,6,81.63,16010,14,914.7,1049,24680,1372,53.3],  // MS derived; check availability
 ["350x350x6.3",66.9,350,6.3,85.25,16640,13.97,951.1,1093,25940,1436,50.6],  // MS derived; check availability
 ["350x350x8.0",84.2,350,8,107,20700,13.9,1180,1370,32600,1790,38.8],
 ["350x350x10.0",104,350,10,133,25200,13.8,1440,1680,40100,2180,30],
 ["350x350x12.0",123,350,12,156,29100,13.6,1660,1950,47600,2550,24.2],
 ["350x350x12.5",127,350,12.5,162,30000,13.6,1720,2020,49400,2640,23],
 ["400x400x6.0",73.3,400,6,93.63,24100,16.04,1205,1379,37040,1808,61.7],  // MS derived; check availability
 ["400x400x6.3",76.8,400,6.3,97.85,25100,16.01,1255,1438,38920,1892,58.5],  // MS derived; check availability
 ["400x400x8.0",96.7,400,8,123,31300,15.9,1560,1800,48900,2360,45],
 ["400x400x10.0",120,400,10,153,38200,15.8,1910,2210,60400,2890,35],
 ["400x400x12.0",141,400,12,180,44300,15.7,2220,2590,71800,3400,28.3],
 ["400x400x12.5",147,400,12.5,187,45900,15.7,2290,2680,74600,3520,27]
].map(a=>({key:a[0],mass:a[1],D:a[2],t:a[3],A:a[4],I:a[5],r:a[6],Z:a[7],S:a[8],J:a[9],C:a[10],dt:a[11]}));

const SHS_HFmap=Object.fromEntries(SHS_HF.map(s=>[s.key,s]));
const SHS_CFmap=Object.fromEntries(SHS_CF.map(s=>[s.key,s]));

const TP385_SHS={"40x40x3.0":[16.5,7.42],"40x40x3.2":[16.5,7.42],"40x40x4.0":[19.5,8.54],"40x40x5.0":[22.5,9.6],"50x50x3.0":[33.8,12.4],"50x50x3.2":[33.8,12.4],"50x50x4.0":[40.4,14.5],"50x50x5.0":[47.6,16.7],"50x50x6.3":[55.2,18.8],"60x60x3.0":[60.2,18.6],"60x60x3.2":[60.2,18.6],"60x60x4.0":[72.5,22],"60x60x5.0":[86.4,25.7],"60x60x6.3":[102,29.6],"60x60x8.0":[118,33.4],"70x70x5.0":[142,36.8],"70x70x6.3":[169,42.9],"70x70x8.0":[200,49.2],"80x80x4.0":[180,41.9],"80x80x5.0":[217,49.8],"80x80x6.3":[262,58.7],"80x80x8.0":[312,68.3],"90x90x4.0":[260,54.2],"90x90x5.0":[316,64.8],"90x90x6.3":[382,77],"90x90x8.0":[459,90.5],"100x100x4.0":[361,68.2],"100x100x5.0":[439,81.8],"100x100x6.3":[534,97.8],"100x100x8.0":[646,116],"100x100x10.0":[761,133],"120x120x5.0":[777,122],"120x120x6.3":[950,147],"120x120x8.0":[1160,176],"120x120x10.0":[1380,206],"120x120x12.5":[1620,236],"140x140x5.0":[1250,170],"140x140x6.3":[1540,206],"140x140x8.0":[1890,249],"140x140x10.0":[2270,294],"140x140x12.5":[2700,342],"150x150x5.0":[1550,197],"150x150x6.3":[1910,240],"150x150x8.0":[2350,291],"150x150x10.0":[2830,344],"150x150x12.5":[3370,402],"160x160x6.3":[2330,275],"160x160x8.0":[2880,335],"160x160x10.0":[3480,398],"160x160x12.5":[4160,467],"180x180x6.3":[3360,355],"180x180x8.0":[4160,434],"180x180x10.0":[5050,518],"180x180x12.5":[6070,613],"180x180x16.0":[7340,724],"200x200x5.0":[3760,362],"200x200x6.3":[4650,444],"200x200x8.0":[5780,545],"200x200x10.0":[7030,655],"200x200x12.5":[8490,778],"200x200x16.0":[10300,927],"250x250x6.3":[9240,712],"250x250x8.0":[11500,880],"250x250x10.0":[14100,1060],"250x250x12.5":[17200,1280],"250x250x16.0":[21100,1550],"300x300x8.0":[20200,1290],"300x300x10.0":[24800,1580],"300x300x12.5":[30300,1900],"300x300x16.0":[37600,2330],"350x350x8.0":[32400,1790],"350x350x10.0":[39900,2190],"350x350x12.5":[48900,2650],"350x350x16.0":[61000,3260],"400x400x10.0":[60100,2900],"400x400x12.5":[73900,3530],"400x400x16.0":[92400,4360],"400x400x20.0":[112000,5240],"70x70x3.0":[98,26.1],"70x70x3.5":[98,26.1],"70x70x4.0":[118,31.2],"80x80x3.0":[148,34.9],"80x80x3.5":[148,34.9],"80x80x6.0":[262,58.7],"90x90x6.0":[382,77],"100x100x6.0":[534,97.8],"120x120x6.0":[950,147],"140x140x6.0":[1540,206],"150x150x6.0":[1910,240],"160x160x6.0":[2330,275],"180x180x6.0":[3360,355],"200x200x6.0":[4650,444],"250x250x6.0":[9240,712]};


/* ==== beam-v03 module 4 ==== */
/* RHS section data and P385 torsion properties. */
const RHS=[
 // key, mass, D(h), B(b), t, A, bT(flange cf/t), dt(web cw/t), Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, J
 // Rows tagged "MS derived": key/D/B/t/mass from the reference software 2025 UK library (EN 10210 hot-finished, S355);
 // A, Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, J derived from the geometry with EN 10210-2:2006 Annex A (ro = 1.5t, ri = 1.0t;
 // J per the ctBoxEN10210() convention, mean corner radius 1.25t). Validated against every section tables row above:
 // max deviation 0.52 % (docs/SECTION_LIBRARY.md). bT = (B - 3t)/t, dt = (D - 3t)/t as this file.
 ["50 x 25 x 2.5",2.69,50.0,25.0,2.5,3.433,7.0,17.0,10.41,3.394,1.741,0.9943,4.163,2.715,5.326,3.222,8.421],  // MS derived
 ["50 x 25 x 3.0",3.17,50.0,25.0,3.0,4.043,5.33,13.7,11.9,3.825,1.716,0.9726,4.762,3.06,6.177,3.71,9.644],  // MS derived
 ["50 x 25 x 3.2",3.36,50.0,25.0,3.2,4.281,4.81,12.6,12.45,3.978,1.706,0.964,4.981,3.182,6.498,3.891,10.09],  // MS derived
 ["50 x 30 x 2.5",2.89,50.0,30.0,2.5,3.683,9.0,17.0,11.82,5.225,1.791,1.191,4.727,3.483,5.92,4.112,11.72],  // MS derived
 ["50 x 30 x 2.9",3.35,50.0,30.0,2.9,4.213,7.34,14.2,13.23,5.804,1.772,1.174,5.292,3.869,6.695,4.633,13.17],  // MS derived
 ["50 x 30 x 3.0",3.41,50.0,30.0,3.0,4.343,7.0,13.7,13.56,5.939,1.767,1.169,5.425,3.959,6.882,4.758,13.52],  // MS derived; check availability
 ["50 x 30 x 3.2",3.61,50.0,30.0,3.2,4.6,6.38,12.6,14.2,6.2,1.76,1.16,5.68,4.13,7.25,5.0,14.2],
 ["50 x 30 x 3.6",4.01,50.0,30.0,3.6,5.1,5.33,10.9,15.4,6.67,1.74,1.14,6.16,4.45,7.94,5.46,15.4],
 ["50 x 30 x 4.0",4.39,50.0,30.0,4.0,5.59,4.5,9.5,16.5,7.08,1.72,1.13,6.6,4.72,8.59,5.88,16.6],
 ["50 x 30 x 5.0",5.28,50.0,30.0,5.0,6.73,3.0,7.0,18.7,7.89,1.67,1.08,7.49,5.26,10.0,6.8,19.0],
 ["50 x 30 x 5.6",5.93,50.0,30.0,5.6,7.369,2.36,5.93,19.73,8.225,1.636,1.056,7.894,5.483,10.76,7.246,20.08],  // MS derived; check availability
 ["50 x 30 x 6.3",6.33,50.0,30.0,6.3,8.07,1.76,4.94,20.6,8.5,1.6,1.03,8.26,5.66,11.5,7.68,21.1],
 ["60 x 40 x 2.5",3.68,60.0,40.0,2.5,4.683,13.0,21.0,22.84,12.06,2.208,1.605,7.612,6.03,9.324,7.016,25.08],  // MS derived
 ["60 x 40 x 2.9",4.26,60.0,40.0,2.9,5.373,10.8,17.7,25.76,13.54,2.189,1.587,8.586,6.769,10.6,7.961,28.42],  // MS derived
 ["60 x 40 x 3.0",4.35,60.0,40.0,3.0,5.543,10.3,17.0,26.46,13.89,2.185,1.583,8.819,6.946,10.91,8.19,29.23],  // MS derived; check availability
 ["60 x 40 x 3.2",4.62,60.0,40.0,3.2,5.88,9.5,15.8,27.8,14.6,2.18,1.57,9.27,7.29,11.5,8.64,30.8],
 ["60 x 40 x 3.6",5.14,60.0,40.0,3.6,6.54,8.11,13.7,30.4,15.9,2.16,1.56,10.1,7.93,12.7,9.5,33.8],
 ["60 x 40 x 4.0",5.64,60.0,40.0,4.0,7.19,7.0,12.0,32.8,17.0,2.14,1.54,10.9,8.52,13.8,10.3,36.7],
 ["60 x 40 x 5.0",6.85,60.0,40.0,5.0,8.73,5.0,9.0,38.1,19.5,2.09,1.5,12.7,9.77,16.4,12.2,43.0],
 ["60 x 40 x 5.6",7.69,60.0,40.0,5.6,9.609,4.14,7.71,40.75,20.75,2.059,1.469,13.58,10.37,17.77,13.14,46.24],  // MS derived; check availability
 ["60 x 40 x 6.0",7.99,60.0,40.0,6.0,10.17,3.67,7.0,42.32,21.45,2.04,1.452,14.11,10.72,18.63,13.73,48.2],  // MS derived
 ["60 x 40 x 6.3",8.31,60.0,40.0,6.3,10.6,3.35,6.52,43.4,21.9,2.02,1.44,14.5,11.0,19.2,14.2,49.5],
 ["60 x 40 x 7.1",9.37,60.0,40.0,7.1,11.64,2.63,5.45,45.85,22.95,1.984,1.404,15.28,11.47,20.7,15.15,52.69],  // MS derived; check availability
 ["60 x 40 x 8.0",10.0,60.0,40.0,8.0,12.8,2.0,4.5,47.9,23.7,1.94,1.36,16.0,11.9,22.1,16.1,55.4],
 ["76.2 x 50.8 x 3.0",5.62,76.2,50.8,3.0,7.163,13.9,22.4,56.67,29.98,2.813,2.046,14.87,11.8,18.17,13.68,62.15],  // MS derived
 ["76.2 x 50.8 x 3.2",5.97,76.2,50.8,3.2,7.609,12.9,20.8,59.79,31.57,2.803,2.037,15.69,12.43,19.23,14.47,65.69],  // MS derived
 ["76.2 x 50.8 x 3.6",6.66,76.2,50.8,3.6,8.487,11.1,18.2,65.8,34.61,2.784,2.02,17.27,13.63,21.3,16.0,72.55],  // MS derived
 ["76.2 x 50.8 x 4.0",7.34,76.2,50.8,4.0,9.348,9.7,16.1,71.5,37.47,2.766,2.002,18.77,14.75,23.29,17.47,79.11],  // MS derived
 ["76.2 x 50.8 x 5.0",8.97,76.2,50.8,5.0,11.43,7.16,12.2,84.45,43.85,2.718,1.959,22.16,17.26,27.97,20.88,94.24],  // MS derived
 ["76.2 x 50.8 x 6.0",10.5,76.2,50.8,6.0,13.41,5.47,9.7,95.6,49.18,2.67,1.915,25.09,19.36,32.2,23.93,107.6],  // MS derived
 ["76.2 x 50.8 x 6.3",11.0,76.2,50.8,6.3,13.99,5.06,9.1,98.61,50.59,2.655,1.902,25.88,19.92,33.39,24.78,111.2],  // MS derived
 ["76.2 x 50.8 x 8.0",13.4,76.2,50.8,8.0,17.07,3.35,6.53,112.9,57.01,2.571,1.827,29.63,22.45,39.39,28.98,128.9],  // MS derived
 ["80 x 40 x 2.9",5.17,80.0,40.0,2.9,6.533,10.8,24.6,52.72,17.54,2.841,1.638,13.18,8.769,16.56,10.11,42.56],  // MS derived
 ["80 x 40 x 3.0",5.29,80.0,40.0,3.0,6.743,10.3,23.7,54.23,18.01,2.836,1.634,13.56,9.004,17.06,10.41,43.79],  // MS derived; check availability
 ["80 x 40 x 3.2",5.62,80.0,40.0,3.2,7.16,9.5,22.0,57.2,18.9,2.83,1.63,14.3,9.46,18.0,11.0,46.2],
 ["80 x 40 x 3.6",6.27,80.0,40.0,3.6,7.98,8.11,19.2,62.8,20.6,2.81,1.61,15.7,10.3,20.0,12.1,50.8],
 ["80 x 40 x 4.0",6.9,80.0,40.0,4.0,8.79,7.0,17.0,68.2,22.2,2.79,1.59,17.1,11.1,21.8,13.2,55.2],
 ["80 x 40 x 5.0",8.42,80.0,40.0,5.0,10.7,5.0,13.0,80.3,25.7,2.74,1.55,20.1,12.9,26.1,15.7,65.1],
 ["80 x 40 x 5.6",9.45,80.0,40.0,5.6,11.85,4.14,11.3,86.65,27.43,2.704,1.522,21.66,13.72,28.5,16.99,70.25],  // MS derived; check availability
 ["80 x 40 x 6.0",9.87,80.0,40.0,6.0,12.57,3.67,10.3,90.54,28.46,2.683,1.504,22.64,14.23,30.0,17.81,73.41],  // MS derived
 ["80 x 40 x 6.3",10.3,80.0,40.0,6.3,13.1,3.35,9.7,93.3,29.2,2.67,1.49,23.3,14.6,31.1,18.4,75.6],
 ["80 x 40 x 7.1",11.6,80.0,40.0,7.1,14.48,2.63,8.27,99.85,30.75,2.626,1.457,24.96,15.37,33.77,19.82,80.91],  // MS derived; check availability
 ["80 x 40 x 8.0",12.5,80.0,40.0,8.0,16.0,2.0,7.0,106.0,32.1,2.58,1.42,26.5,16.1,36.5,21.2,85.8],
 ["90 x 50 x 3.0",6.24,90.0,50.0,3.0,7.943,13.7,27.0,84.38,33.47,3.259,2.053,18.75,13.39,23.19,15.34,76.49],  // MS derived; check availability
 ["90 x 50 x 3.2",6.63,90.0,50.0,3.2,8.44,12.6,25.1,89.1,35.3,3.25,2.04,19.8,14.1,24.6,16.2,80.9],
 ["90 x 50 x 3.6",7.4,90.0,50.0,3.6,9.42,10.9,22.0,98.3,38.7,3.23,2.03,21.8,15.5,27.2,18.0,89.4],
 ["90 x 50 x 4.0",8.15,90.0,50.0,4.0,10.4,9.5,19.5,107.0,41.9,3.21,2.01,23.8,16.8,29.8,19.6,97.5],
 ["90 x 50 x 5.0",9.99,90.0,50.0,5.0,12.7,7.0,15.0,127.0,49.2,3.16,1.97,28.3,19.7,36.0,23.5,116.0],
 ["90 x 50 x 5.6",11.2,90.0,50.0,5.6,14.09,5.93,13.1,138.2,53.03,3.132,1.94,30.71,21.21,39.43,25.68,126.7],  // MS derived; check availability
 ["90 x 50 x 6.0",11.8,90.0,50.0,6.0,14.97,5.33,12.0,145.0,55.36,3.112,1.923,32.22,22.14,41.63,27.04,133.1],  // MS derived
 ["90 x 50 x 6.3",12.3,90.0,50.0,6.3,15.6,4.94,11.3,150.0,57.0,3.1,1.91,33.3,22.8,43.2,28.0,138.0],
 ["90 x 50 x 7.1",13.8,90.0,50.0,7.1,17.32,4.04,9.68,161.8,60.91,3.056,1.875,35.96,24.36,47.25,30.47,149.1],  // MS derived; check availability
 ["90 x 50 x 8.0",15.0,90.0,50.0,8.0,19.2,3.25,8.25,174.0,64.6,3.01,1.84,38.6,25.8,51.4,32.9,160.0],
 ["90 x 50 x 8.8",16.6,90.0,50.0,8.8,20.71,2.68,7.23,182.6,67.23,2.969,1.802,40.57,26.89,54.78,34.9,168.9],  // MS derived; check availability
 ["90 x 50 x 10.0",18.0,90.0,50.0,10.0,22.9,2.0,6.0,194.0,70.2,2.91,1.75,43.0,28.1,59.3,37.4,179.0],
 ["100 x 50 x 3.0",6.71,100.0,50.0,3.0,8.543,13.7,30.3,109.6,36.79,3.582,2.075,21.92,14.72,27.31,16.75,88.37],  // MS derived
 ["100 x 50 x 3.2",7.13,100.0,50.0,3.2,9.08,12.6,28.3,116.0,38.8,3.57,2.07,23.2,15.5,28.9,17.7,93.4],
 ["100 x 50 x 3.6",7.96,100.0,50.0,3.6,10.1,10.9,24.8,128.0,42.6,3.55,2.05,25.6,17.0,32.1,19.6,103.0],
 ["100 x 50 x 4.0",8.78,100.0,50.0,4.0,11.2,9.5,22.0,140.0,46.2,3.53,2.03,27.9,18.5,35.2,21.5,113.0],
 ["100 x 50 x 5.0",10.8,100.0,50.0,5.0,13.7,7.0,17.0,167.0,54.3,3.48,1.99,33.3,21.7,42.6,25.8,135.0],
 ["100 x 50 x 5.6",12.1,100.0,50.0,5.6,15.21,5.93,14.9,181.2,58.58,3.452,1.963,36.25,23.43,46.76,28.17,146.8],  // MS derived; check availability
 ["100 x 50 x 6.0",12.7,100.0,50.0,6.0,16.17,5.33,13.7,190.5,61.2,3.432,1.945,38.09,24.48,49.41,29.68,154.3],  // MS derived
 ["100 x 50 x 6.3",13.3,100.0,50.0,6.3,16.9,4.94,12.9,197.0,63.0,3.42,1.93,39.4,25.2,51.3,30.8,160.0],
 ["100 x 50 x 7.1",14.9,100.0,50.0,7.1,18.74,4.04,11.1,213.5,67.5,3.375,1.898,42.7,27.0,56.26,33.51,173.1],  // MS derived; check availability
 ["100 x 50 x 8.0",16.3,100.0,50.0,8.0,20.8,3.25,9.5,230.0,71.7,3.33,1.86,46.0,28.7,61.4,36.3,186.0],
 ["100 x 50 x 8.8",18.0,100.0,50.0,8.8,22.47,2.68,8.36,242.7,74.81,3.286,1.825,48.53,29.92,65.57,38.52,196.7],  // MS derived; check availability
 ["100 x 50 x 10.0",19.6,100.0,50.0,10.0,24.9,2.0,7.0,259.0,78.4,3.22,1.77,51.8,31.4,71.2,41.4,209.0],
 ["100 x 60 x 3.0",7.18,100.0,60.0,3.0,9.143,17.0,30.3,123.7,55.73,3.678,2.469,24.74,18.58,30.22,21.17,121.4],  // MS derived
 ["100 x 60 x 3.2",7.63,100.0,60.0,3.2,9.72,15.8,28.3,131.0,58.8,3.67,2.46,26.2,19.6,32.0,22.4,129.0],
 ["100 x 60 x 3.6",8.53,100.0,60.0,3.6,10.9,13.7,24.8,145.0,64.8,3.65,2.44,28.9,21.6,35.6,24.9,142.0],
 ["100 x 60 x 4.0",9.41,100.0,60.0,4.0,12.0,12.0,22.0,158.0,70.5,3.63,2.43,31.6,23.5,39.1,27.3,156.0],
 ["100 x 60 x 5.0",11.6,100.0,60.0,5.0,14.7,9.0,17.0,189.0,83.6,3.58,2.38,37.8,27.9,47.4,32.9,188.0],
 ["100 x 60 x 5.6",13.0,100.0,60.0,5.6,16.33,7.71,14.9,206.2,90.64,3.554,2.356,41.24,30.21,52.04,36.05,205.1],  // MS derived; check availability
 ["100 x 60 x 6.0",13.6,100.0,60.0,6.0,17.37,7.0,13.7,217.0,95.02,3.534,2.339,43.4,31.67,55.05,38.07,216.3],  // MS derived
 ["100 x 60 x 6.3",14.2,100.0,60.0,6.3,18.1,6.52,12.9,225.0,98.1,3.52,2.33,45.0,32.7,57.3,39.5,224.0],
 ["100 x 60 x 7.1",16.1,100.0,60.0,7.1,20.16,5.45,11.1,244.2,105.8,3.48,2.291,48.84,35.27,62.86,43.24,244.7],  // MS derived; check availability
 ["100 x 60 x 8.0",17.5,100.0,60.0,8.0,22.4,4.5,9.5,264.0,113.0,3.44,2.25,52.8,37.8,68.7,47.1,265.0],
 ["100 x 60 x 8.8",19.4,100.0,60.0,8.8,24.23,3.82,8.36,279.4,119.1,3.395,2.217,55.87,39.7,73.6,50.2,282.0],  // MS derived; check availability
 ["100 x 60 x 10.0",21.1,100.0,60.0,10.0,26.9,3.0,7.0,299.0,126.0,3.33,2.16,59.9,42.1,80.2,54.4,304.0],
 ["120 x 60 x 3.6",9.66,120.0,60.0,3.6,12.3,13.7,30.3,227.0,76.3,4.3,2.49,37.9,25.4,47.2,28.9,183.0],
 ["120 x 60 x 4.0",10.7,120.0,60.0,4.0,13.6,12.0,27.0,249.0,83.1,4.28,2.47,41.5,27.7,51.9,31.7,201.0],
 ["120 x 60 x 5.0",13.1,120.0,60.0,5.0,16.7,9.0,21.0,299.0,98.8,4.23,2.43,49.9,32.9,63.1,38.4,242.0],
 ["120 x 60 x 5.6",14.7,120.0,60.0,5.6,18.57,7.71,18.4,327.4,107.3,4.199,2.404,54.56,35.76,69.49,42.14,264.8],  // MS derived; check availability
 ["120 x 60 x 6.0",15.5,120.0,60.0,6.0,19.77,7.0,17.0,345.3,112.6,4.179,2.386,57.55,37.53,73.63,44.55,279.4],  // MS derived
 ["120 x 60 x 6.3",16.2,120.0,60.0,6.3,20.7,6.52,16.0,358.0,116.0,4.16,2.37,59.7,38.8,76.7,46.3,290.0],
 ["120 x 60 x 7.1",18.3,120.0,60.0,7.1,23.0,5.45,13.9,391.0,125.8,4.123,2.339,65.17,41.93,84.44,50.75,316.8],  // MS derived
 ["120 x 60 x 8.0",20.1,120.0,60.0,8.0,25.6,4.5,12.0,425.0,135.0,4.08,2.3,70.8,45.0,92.7,55.4,344.0],
 ["120 x 60 x 8.8",22.1,120.0,60.0,8.8,27.75,3.82,10.6,452.0,142.4,4.036,2.265,75.33,47.47,99.59,59.21,366.4],  // MS derived; check availability
 ["120 x 60 x 10.0",24.3,120.0,60.0,10.0,30.9,3.0,9.0,488.0,152.0,3.97,2.21,81.4,50.5,109.0,64.4,396.0],
 ["120 x 80 x 3.6",10.8,120.0,80.0,3.6,13.7,19.2,30.3,276.0,147.0,4.48,3.27,46.0,36.7,55.6,42.0,301.0],
 ["120 x 80 x 4.0",11.9,120.0,80.0,4.0,15.2,17.0,27.0,303.0,161.0,4.46,3.25,50.4,40.2,61.2,46.1,330.0],
 ["120 x 80 x 5.0",14.7,120.0,80.0,5.0,18.7,13.0,21.0,365.0,193.0,4.42,3.21,60.9,48.2,74.6,56.1,401.0],
 ["120 x 80 x 5.6",16.5,120.0,80.0,5.6,20.81,11.3,18.4,400.7,210.9,4.388,3.183,66.79,52.72,82.31,61.83,441.6],  // MS derived; check availability
 ["120 x 80 x 6.0",17.4,120.0,80.0,6.0,22.17,10.3,17.0,423.3,222.3,4.369,3.166,70.56,55.56,87.31,65.52,467.6],  // MS derived
 ["120 x 80 x 6.3",18.2,120.0,80.0,6.3,23.2,9.7,16.0,440.0,230.0,4.36,3.15,73.3,57.6,91.0,68.2,487.0],
 ["120 x 80 x 7.1",20.5,120.0,80.0,7.1,25.84,8.27,13.9,481.7,251.3,4.317,3.118,80.28,62.81,100.5,75.17,535.3],  // MS derived
 ["120 x 80 x 8.0",22.6,120.0,80.0,8.0,28.8,7.0,12.0,525.0,273.0,4.27,3.08,87.5,68.1,111.0,82.6,587.0],
 ["120 x 80 x 8.8",24.9,120.0,80.0,8.8,31.27,6.09,10.6,561.0,289.7,4.236,3.044,93.5,72.43,119.2,88.72,629.1],  // MS derived; check availability
 ["120 x 80 x 10.0",27.4,120.0,80.0,10.0,34.9,5.0,9.0,609.0,313.0,4.18,2.99,102.0,78.1,131.0,97.3,688.0],
 ["150 x 100 x 4.0",15.1,150.0,100.0,4.0,19.2,22.0,34.5,607.0,324.0,5.63,4.11,81.0,64.8,97.4,73.6,660.0],
 ["150 x 100 x 5.0",18.6,150.0,100.0,5.0,23.7,17.0,27.0,739.0,392.0,5.58,4.07,98.5,78.5,119.0,90.1,807.0],
 ["150 x 100 x 5.6",20.9,150.0,100.0,5.6,26.41,14.9,23.8,813.8,431.0,5.551,4.04,108.5,86.21,132.2,99.62,891.2],  // MS derived; check availability
 ["150 x 100 x 6.0",22.1,150.0,100.0,6.0,28.17,13.7,22.0,862.3,455.9,5.532,4.023,115.0,91.18,140.5,105.8,946.1],  // MS derived
 ["150 x 100 x 6.3",23.1,150.0,100.0,6.3,29.5,12.9,20.8,898.0,474.0,5.52,4.01,120.0,94.8,147.0,110.0,986.0],
 ["150 x 100 x 7.1",26.1,150.0,100.0,7.1,32.94,11.1,18.1,989.5,520.5,5.481,3.975,131.9,104.1,162.7,122.2,1091.0],  // MS derived; check availability
 ["150 x 100 x 8.0",28.9,150.0,100.0,8.0,36.8,9.5,15.8,1090.0,569.0,5.44,3.94,145.0,114.0,180.0,135.0,1200.0],
 ["150 x 100 x 8.8",31.8,150.0,100.0,8.8,40.07,8.36,14.0,1168.0,609.8,5.4,3.901,155.8,122.0,194.9,145.8,1298.0],  // MS derived; check availability
 ["150 x 100 x 10.0",35.3,150.0,100.0,10.0,44.9,7.0,12.0,1280.0,665.0,5.34,3.85,171.0,133.0,216.0,161.0,1430.0],
 ["150 x 100 x 11.0",38.9,150.0,100.0,11.0,48.86,6.09,10.6,1370.0,707.4,5.294,3.805,182.6,141.5,232.7,173.3,1536.0],  // MS derived; check availability
 ["150 x 100 x 12.0",41.4,150.0,100.0,12.0,52.69,5.33,9.5,1450.0,745.4,5.246,3.761,193.4,149.1,248.6,184.6,1633.0],  // MS derived
 ["150 x 100 x 12.5",42.8,150.0,100.0,12.5,54.6,5.0,9.0,1490.0,763.0,5.22,3.74,198.0,153.0,256.0,190.0,1680.0],
 ["160 x 80 x 4.0",14.4,160.0,80.0,4.0,18.4,17.0,37.0,612.0,207.0,5.77,3.35,76.5,51.7,94.7,58.3,493.0],
 ["160 x 80 x 5.0",17.8,160.0,80.0,5.0,22.7,13.0,29.0,744.0,249.0,5.72,3.31,93.0,62.3,116.0,71.1,600.0],
 ["160 x 80 x 5.6",20.0,160.0,80.0,5.6,25.29,11.3,25.6,819.2,273.0,5.691,3.286,102.4,68.25,128.4,78.5,661.1],  // MS derived; check availability
 ["160 x 80 x 6.0",21.2,160.0,80.0,6.0,26.97,10.3,23.7,867.7,288.1,5.672,3.268,108.5,72.03,136.5,83.28,700.6],  // MS derived
 ["160 x 80 x 6.3",22.2,160.0,80.0,6.3,28.2,9.7,22.4,903.0,299.0,5.66,3.26,113.0,74.8,142.0,86.8,730.0],
 ["160 x 80 x 7.1",25.0,160.0,80.0,7.1,31.52,8.27,19.5,994.5,327.0,5.617,3.221,124.3,81.74,157.8,95.88,804.0],  // MS derived; check availability
 ["160 x 80 x 8.0",27.6,160.0,80.0,8.0,35.2,7.0,17.0,1090.0,356.0,5.57,3.18,136.0,89.0,175.0,106.0,883.0],
 ["160 x 80 x 8.8",30.4,160.0,80.0,8.8,38.31,6.09,15.2,1172.0,379.4,5.531,3.147,146.5,94.85,188.7,113.8,949.1],  // MS derived; check availability
 ["160 x 80 x 10.0",33.7,160.0,80.0,10.0,42.9,5.0,13.0,1280.0,411.0,5.47,3.1,161.0,103.0,209.0,125.0,1040.0],
 ["160 x 80 x 11.0",37.2,160.0,80.0,11.0,46.66,4.27,11.5,1370.0,434.6,5.419,3.052,171.3,108.7,224.9,134.2,1111.0],  // MS derived; check availability
 ["160 x 80 x 12.0",39.5,160.0,80.0,12.0,50.29,3.67,10.3,1449.0,455.3,5.367,3.009,181.1,113.8,240.0,142.5,1175.0],  // MS derived
 ["160 x 80 x 12.5",40.9,160.0,80.0,12.5,52.1,3.4,9.8,1480.0,465.0,5.34,2.99,186.0,116.0,247.0,146.0,1200.0],
 ["180 x 60 x 4.0",14.5,180.0,60.0,4.0,18.39,12.0,42.0,696.7,120.8,6.155,2.563,77.41,40.26,99.84,45.19,340.8],  // MS derived; check availability
 ["180 x 60 x 5.0",18.0,180.0,60.0,5.0,22.73,9.0,33.0,846.3,144.3,6.102,2.519,94.04,48.09,122.3,54.89,411.4],  // MS derived; check availability
 ["180 x 60 x 5.6",20.0,180.0,60.0,5.6,25.29,7.71,29.1,931.6,157.2,6.07,2.493,103.5,52.39,135.3,60.42,451.1],  // MS derived; check availability
 ["180 x 60 x 6.3",22.3,180.0,60.0,6.3,28.23,6.52,25.6,1027.0,171.1,6.032,2.462,114.1,57.05,150.0,66.59,494.9],  // MS derived; check availability
 ["180 x 60 x 7.1",25.0,180.0,60.0,7.1,31.52,5.45,22.4,1130.0,185.8,5.988,2.428,125.6,61.92,166.2,73.29,541.7],  // MS derived; check availability
 ["180 x 60 x 8.0",27.9,180.0,60.0,8.0,35.15,4.5,19.5,1240.0,200.6,5.938,2.389,137.7,66.85,183.8,80.36,590.2],  // MS derived; check availability
 ["180 x 60 x 8.8",30.4,180.0,60.0,8.8,38.31,3.82,17.5,1331.0,212.3,5.894,2.354,147.9,70.76,198.7,86.24,629.7],  // MS derived; check availability
 ["180 x 60 x 10.0",34.2,180.0,60.0,10.0,42.93,3.0,15.0,1457.0,227.5,5.827,2.302,161.9,75.85,219.9,94.38,682.8],  // MS derived; check availability
 ["180 x 60 x 11.0",37.2,180.0,60.0,11.0,46.66,2.45,13.4,1554.0,238.2,5.77,2.26,172.6,79.41,236.6,100.5,721.6],  // MS derived; check availability
 ["180 x 60 x 12.5",41.6,180.0,60.0,12.5,52.07,1.8,11.4,1682.0,251.1,5.684,2.196,186.9,83.71,259.9,108.7,770.5],  // MS derived; check availability
 ["180 x 100 x 4.0",17.0,180.0,100.0,4.0,21.59,22.0,42.0,944.5,379.4,6.614,4.192,104.9,75.87,128.0,85.16,852.1],  // MS derived; check availability
 ["180 x 100 x 5.0",21.1,180.0,100.0,5.0,26.73,17.0,33.0,1153.0,460.1,6.567,4.149,128.1,92.02,157.3,104.4,1042.0],  // MS derived; check availability
 ["180 x 100 x 5.6",23.5,180.0,100.0,5.6,29.77,14.9,29.1,1272.0,506.0,6.538,4.123,141.4,101.2,174.3,115.5,1152.0],  // MS derived; check availability
 ["180 x 100 x 6.3",26.3,180.0,100.0,6.3,33.27,12.9,25.6,1407.0,557.2,6.504,4.092,156.4,111.4,193.8,128.1,1277.0],  // MS derived
 ["180 x 100 x 7.1",29.4,180.0,100.0,7.1,37.2,11.1,22.4,1555.0,612.6,6.465,4.058,172.8,122.5,215.3,142.0,1413.0],  // MS derived; check availability
 ["180 x 100 x 8.0",32.9,180.0,100.0,8.0,41.55,9.5,19.5,1713.0,671.1,6.421,4.019,190.4,134.2,238.8,157.1,1560.0],  // MS derived
 ["180 x 100 x 8.8",36.0,180.0,100.0,8.8,45.35,8.36,17.5,1847.0,719.9,6.382,3.984,205.3,144.0,258.9,169.9,1685.0],  // MS derived; check availability
 ["180 x 100 x 10.0",40.4,180.0,100.0,10.0,50.93,7.0,15.0,2036.0,787.4,6.323,3.932,226.2,157.5,287.9,188.2,1862.0],  // MS derived
 ["180 x 100 x 11.0",44.1,180.0,100.0,11.0,55.46,6.09,13.4,2183.0,838.7,6.273,3.889,242.5,167.7,311.0,202.7,2000.0],  // MS derived; check availability
 ["180 x 100 x 12.5",49.5,180.0,100.0,12.5,62.07,5.0,11.4,2385.0,907.6,6.198,3.824,265.0,181.5,343.7,222.9,2191.0],  // MS derived
 ["180 x 100 x 14.2",55.3,180.0,100.0,14.2,69.29,4.04,9.68,2589.0,974.5,6.113,3.75,287.7,194.9,378.0,243.7,2385.0],  // MS derived; check availability
 ["200 x 100 x 4.0",18.2,200.0,100.0,4.0,23.2,22.0,47.0,1220.0,416.0,7.26,4.24,122.0,83.2,150.0,92.8,983.0],
 ["200 x 100 x 5.0",22.6,200.0,100.0,5.0,28.7,17.0,37.0,1500.0,505.0,7.21,4.19,149.0,101.0,185.0,114.0,1200.0],
 ["200 x 100 x 5.6",25.3,200.0,100.0,5.6,32.01,14.9,32.7,1652.0,555.9,7.183,4.168,165.2,111.2,205.2,126.1,1331.0],  // MS derived; check availability
 ["200 x 100 x 6.0",26.8,200.0,100.0,6.0,34.17,13.7,30.3,1754.0,588.6,7.164,4.15,175.4,117.7,218.5,134.0,1414.0],  // MS derived
 ["200 x 100 x 6.3",28.1,200.0,100.0,6.3,35.8,12.9,28.7,1830.0,613.0,7.15,4.14,183.0,123.0,228.0,140.0,1480.0],
 ["200 x 100 x 7.1",31.7,200.0,100.0,7.1,40.04,11.1,25.2,2024.0,674.0,7.109,4.103,202.4,134.8,254.0,155.2,1634.0],  // MS derived; check availability
 ["200 x 100 x 8.0",35.1,200.0,100.0,8.0,44.8,9.5,22.0,2230.0,739.0,7.06,4.06,223.0,148.0,282.0,172.0,1800.0],
 ["200 x 100 x 8.8",38.7,200.0,100.0,8.8,48.87,8.36,19.7,2412.0,793.3,7.025,4.029,241.2,158.7,306.1,186.0,1950.0],  // MS derived; check availability
 ["200 x 100 x 10.0",43.1,200.0,100.0,10.0,54.9,7.0,17.0,2660.0,869.0,6.96,3.98,266.0,174.0,341.0,206.0,2160.0],
 ["200 x 100 x 11.0",47.6,200.0,100.0,11.0,59.86,6.09,15.2,2862.0,926.3,6.914,3.934,286.2,185.3,368.6,222.2,2317.0],  // MS derived; check availability
 ["200 x 100 x 12.0",50.8,200.0,100.0,12.0,64.69,5.33,13.7,3047.0,979.2,6.863,3.89,304.7,195.8,395.3,237.4,2469.0],  // MS derived
 ["200 x 100 x 12.5",52.7,200.0,100.0,12.5,67.1,5.0,13.0,3140.0,1000.0,6.84,3.87,314.0,201.0,408.0,245.0,2540.0],
 ["200 x 100 x 14.2",58.9,200.0,100.0,14.2,75.0,4.04,11.1,3420.0,1080.0,6.75,3.8,342.0,216.0,450.0,268.0,2770.0],
 ["200 x 100 x 16.0",65.2,200.0,100.0,16.0,83.0,3.25,9.5,3680.0,1150.0,6.66,3.72,368.0,229.0,491.0,290.0,2980.0],
 ["200 x 120 x 5.0",24.1,200.0,120.0,5.0,30.7,21.0,37.0,1680.0,762.0,7.4,4.98,168.0,127.0,205.0,144.0,1650.0],
 ["200 x 120 x 5.6",27.0,200.0,120.0,5.6,34.25,18.4,32.7,1863.0,840.8,7.376,4.955,186.3,140.1,227.0,159.2,1826.0],  // MS derived; check availability
 ["200 x 120 x 6.0",28.7,200.0,120.0,6.0,36.57,17.0,30.3,1980.0,891.6,7.357,4.937,198.0,148.6,241.8,169.4,1942.0],  // MS derived
 ["200 x 120 x 6.3",30.1,200.0,120.0,6.3,38.3,16.0,28.7,2060.0,929.0,7.34,4.92,207.0,155.0,253.0,177.0,2030.0],
 ["200 x 120 x 7.1",33.9,200.0,120.0,7.1,42.88,13.9,25.2,2288.0,1025.0,7.305,4.89,228.8,170.9,281.3,196.7,2252.0],  // MS derived; check availability
 ["200 x 120 x 8.0",37.6,200.0,120.0,8.0,48.0,12.0,22.0,2530.0,1130.0,7.26,4.85,253.0,188.0,313.0,218.0,2500.0],
 ["200 x 120 x 8.8",41.5,200.0,120.0,8.8,52.39,10.6,19.7,2734.0,1215.0,7.223,4.816,273.4,202.5,339.7,236.6,2703.0],  // MS derived; check availability
 ["200 x 120 x 10.0",46.3,200.0,120.0,10.0,58.9,9.0,17.0,3030.0,1340.0,7.17,4.76,303.0,223.0,379.0,263.0,3000.0],
 ["200 x 120 x 11.0",51.0,200.0,120.0,11.0,64.26,7.91,15.2,3255.0,1432.0,7.117,4.721,325.5,238.7,410.2,284.3,3237.0],  // MS derived; check availability
 ["200 x 120 x 12.5",56.6,200.0,120.0,12.5,72.1,6.6,13.0,3580.0,1560.0,7.04,4.66,358.0,260.0,455.0,314.0,3570.0],
 ["200 x 120 x 14.2",63.3,200.0,120.0,14.2,80.7,5.45,11.1,3910.0,1690.0,6.96,4.58,391.0,282.0,503.0,346.0,3920.0],
 ["200 x 120 x 16.0",70.2,200.0,120.0,16.0,89.4,4.5,9.5,4220.0,1810.0,6.87,4.5,422.0,302.0,550.0,377.0,4250.0],
 ["200 x 150 x 5.0",26.5,200.0,150.0,5.0,33.7,27.0,37.0,1970.0,1260.0,7.64,6.12,197.0,169.0,234.0,192.0,2390.0],
 ["200 x 150 x 5.6",29.7,200.0,150.0,5.6,37.61,23.8,32.7,2181.0,1398.0,7.615,6.097,218.1,186.4,259.7,213.1,2648.0],  // MS derived; check availability
 ["200 x 150 x 6.3",33.0,200.0,150.0,6.3,42.1,20.8,28.7,2420.0,1550.0,7.58,6.07,242.0,207.0,289.0,237.0,2950.0],
 ["200 x 150 x 7.1",37.2,200.0,150.0,7.1,47.14,18.1,25.2,2685.0,1715.0,7.546,6.032,268.5,228.7,322.4,264.2,3280.0],  // MS derived; check availability
 ["200 x 150 x 8.0",41.4,200.0,150.0,8.0,52.8,15.8,22.0,2970.0,1890.0,7.5,5.99,297.0,253.0,359.0,294.0,3640.0],
 ["200 x 150 x 8.8",45.6,200.0,150.0,8.8,57.67,14.0,19.7,3217.0,2047.0,7.468,5.958,321.7,272.9,390.2,319.1,3956.0],  // MS derived; check availability
 ["200 x 150 x 10.0",51.0,200.0,150.0,10.0,64.9,12.0,17.0,3570.0,2260.0,7.41,5.91,357.0,302.0,436.0,356.0,4410.0],
 ["200 x 150 x 11.0",56.2,200.0,150.0,11.0,70.86,10.6,15.2,3845.0,2435.0,7.366,5.861,384.5,324.6,472.6,385.6,4771.0],  // MS derived; check availability
 ["200 x 150 x 12.5",62.5,200.0,150.0,12.5,79.6,9.0,13.0,4240.0,2670.0,7.3,5.8,424.0,356.0,525.0,428.0,5290.0],
 ["200 x 150 x 14.2",70.0,200.0,150.0,14.2,89.2,7.56,11.1,4640.0,2920.0,7.22,5.72,464.0,389.0,582.0,473.0,5830.0],
 ["200 x 150 x 16.0",77.7,200.0,150.0,16.0,99.0,6.38,9.5,5040.0,3150.0,7.13,5.64,504.0,420.0,638.0,518.0,6370.0],
 ["220 x 120 x 7.1",36.1,220.0,120.0,7.1,45.72,13.9,28.0,2895.0,1116.0,7.957,4.94,263.2,186.0,325.7,212.7,2572.0],  // MS derived; check availability
 ["220 x 120 x 8.0",40.5,220.0,120.0,8.0,51.15,12.0,24.5,3203.0,1229.0,7.913,4.901,291.2,204.8,362.2,236.1,2850.0],  // MS derived
 ["220 x 120 x 8.8",44.2,220.0,120.0,8.8,55.91,10.6,22.0,3467.0,1324.0,7.874,4.867,315.1,220.7,393.9,256.2,3089.0],  // MS derived; check availability
 ["220 x 120 x 10.0",49.9,220.0,120.0,10.0,62.93,9.0,19.0,3844.0,1459.0,7.815,4.815,349.4,243.1,439.8,285.2,3431.0],  // MS derived
 ["220 x 120 x 11.0",54.5,220.0,120.0,11.0,68.66,7.91,17.0,4141.0,1563.0,7.766,4.772,376.5,260.5,476.7,308.3,3703.0],  // MS derived; check availability
 ["220 x 120 x 12.5",61.2,220.0,120.0,12.5,77.07,6.6,14.6,4560.0,1707.0,7.692,4.707,414.5,284.5,529.7,341.2,4087.0],  // MS derived
 ["220 x 120 x 14.2",68.7,220.0,120.0,14.2,86.33,5.45,12.5,4996.0,1853.0,7.607,4.633,454.1,308.8,586.4,376.0,4488.0],  // MS derived; check availability
 ["220 x 120 x 16.0",76.4,220.0,120.0,16.0,95.81,4.5,10.8,5413.0,1988.0,7.516,4.555,492.1,331.3,642.6,409.9,4873.0],  // MS derived
 ["250 x 100 x 5.0",26.5,250.0,100.0,5.0,33.7,17.0,47.0,2610.0,618.0,8.8,4.28,209.0,124.0,263.0,138.0,1620.0],
 ["250 x 100 x 5.6",29.7,250.0,100.0,5.6,37.61,14.9,41.6,2889.0,680.9,8.765,4.255,231.2,136.2,292.3,152.5,1789.0],  // MS derived; check availability
 ["250 x 100 x 6.3",33.0,250.0,100.0,6.3,42.1,12.9,36.7,3210.0,751.0,8.73,4.22,257.0,150.0,326.0,169.0,1980.0],
 ["250 x 100 x 7.1",37.2,250.0,100.0,7.1,47.14,11.1,32.2,3559.0,827.5,8.688,4.19,284.7,165.5,362.9,188.2,2198.0],  // MS derived; check availability
 ["250 x 100 x 8.0",41.4,250.0,100.0,8.0,52.8,9.5,28.3,3940.0,909.0,8.64,4.15,315.0,182.0,404.0,209.0,2430.0],
 ["250 x 100 x 8.8",45.6,250.0,100.0,8.8,57.67,8.36,25.4,4266.0,976.9,8.6,4.116,341.3,195.4,439.2,226.1,2627.0],  // MS derived; check availability
 ["250 x 100 x 10.0",51.0,250.0,100.0,10.0,64.9,7.0,22.0,4730.0,1070.0,8.54,4.06,379.0,214.0,491.0,251.0,2910.0],
 ["250 x 100 x 11.0",56.2,250.0,100.0,11.0,70.86,6.09,19.7,5102.0,1145.0,8.485,4.02,408.2,229.0,532.0,271.2,3128.0],  // MS derived; check availability
 ["250 x 100 x 12.5",62.5,250.0,100.0,12.5,79.6,5.0,17.0,5620.0,1240.0,8.41,3.96,450.0,249.0,592.0,299.0,3440.0],
 ["250 x 100 x 14.2",70.0,250.0,100.0,14.2,89.2,4.04,14.6,6160.0,1340.0,8.31,3.88,493.0,269.0,655.0,329.0,3750.0],
 ["250 x 100 x 16.0",77.7,250.0,100.0,16.0,99.0,3.25,12.6,6690.0,1430.0,8.22,3.8,535.0,287.0,719.0,358.0,4050.0],
 ["250 x 150 x 5.0",30.4,250.0,150.0,5.0,38.7,27.0,47.0,3360.0,1530.0,9.31,6.28,269.0,204.0,324.0,228.0,3280.0],
 ["250 x 150 x 5.6",34.1,250.0,150.0,5.6,43.21,23.8,41.6,3726.0,1690.0,9.286,6.254,298.1,225.3,360.7,253.5,3640.0],  // MS derived; check availability
 ["250 x 150 x 6.0",36.2,250.0,150.0,6.0,46.17,22.0,38.7,3965.0,1796.0,9.267,6.237,317.2,239.5,384.6,270.1,3877.0],  // MS derived
 ["250 x 150 x 6.3",38.0,250.0,150.0,6.3,48.4,20.8,36.7,4140.0,1870.0,9.25,6.22,331.0,250.0,402.0,283.0,4050.0],
 ["250 x 150 x 7.1",42.8,250.0,150.0,7.1,54.24,18.1,32.2,4606.0,2078.0,9.215,6.189,368.5,277.0,449.2,314.9,4515.0],  // MS derived; check availability
 ["250 x 150 x 8.0",47.7,250.0,150.0,8.0,60.8,15.8,28.3,5110.0,2300.0,9.17,6.15,409.0,306.0,501.0,350.0,5020.0],
 ["250 x 150 x 8.8",52.5,250.0,150.0,8.8,66.47,14.0,25.4,5546.0,2486.0,9.134,6.116,443.7,331.5,545.4,381.3,5457.0],  // MS derived; check availability
 ["250 x 150 x 10.0",58.8,250.0,150.0,10.0,74.9,12.0,22.0,6170.0,2760.0,9.08,6.06,494.0,367.0,611.0,426.0,6090.0],
 ["250 x 150 x 11.0",64.8,250.0,150.0,11.0,81.86,10.6,19.7,6674.0,2967.0,9.029,6.02,533.9,395.6,663.5,462.1,6598.0],  // MS derived; check availability
 ["250 x 150 x 12.0",69.6,250.0,150.0,12.0,88.69,9.5,17.8,7154.0,3168.0,8.981,5.977,572.3,422.5,714.8,497.0,7088.0],  // MS derived
 ["250 x 150 x 12.5",72.3,250.0,150.0,12.5,92.1,9.0,17.0,7390.0,3260.0,8.96,5.96,591.0,435.0,740.0,514.0,7330.0],
 ["250 x 150 x 14.2",81.1,250.0,150.0,14.2,103.0,7.56,14.6,8140.0,3580.0,8.87,5.88,651.0,477.0,823.0,570.0,8100.0],
 ["250 x 150 x 16.0",90.3,250.0,150.0,16.0,115.0,6.38,12.6,8880.0,3870.0,8.79,5.8,710.0,516.0,906.0,625.0,8870.0],
 ["250 x 150 x 17.5",97.7,250.0,150.0,17.5,124.0,5.57,11.3,9450.0,4100.0,8.71,5.74,756.0,546.0,972.0,669.0,9460.0],
 ["260 x 140 x 5.0",30.5,260.0,140.0,5.0,38.73,25.0,49.0,3532.0,1354.0,9.549,5.914,271.7,193.5,331.2,215.8,3078.0],  // MS derived
 ["260 x 140 x 5.6",34.1,260.0,140.0,5.6,43.21,22.0,43.4,3916.0,1498.0,9.52,5.888,301.2,214.0,368.3,239.7,3417.0],  // MS derived; check availability
 ["260 x 140 x 6.3",38.0,260.0,140.0,6.3,48.39,19.2,38.3,4355.0,1660.0,9.487,5.857,335.0,237.2,410.9,267.0,3803.0],  // MS derived
 ["260 x 140 x 7.1",42.8,260.0,140.0,7.1,54.24,16.7,33.6,4842.0,1839.0,9.448,5.823,372.5,262.7,458.7,297.6,4234.0],  // MS derived; check availability
 ["260 x 140 x 8.0",47.7,260.0,140.0,8.0,60.75,14.5,29.5,5373.0,2032.0,9.404,5.784,413.3,290.3,511.3,331.1,4704.0],  // MS derived
 ["260 x 140 x 8.8",52.5,260.0,140.0,8.8,66.47,12.9,26.5,5831.0,2197.0,9.366,5.749,448.5,313.9,556.9,360.0,5110.0],  // MS derived; check availability
 ["260 x 140 x 10.0",58.8,260.0,140.0,10.0,74.93,11.0,23.0,6490.0,2432.0,9.307,5.697,499.3,347.4,623.6,402.1,5698.0],  // MS derived
 ["260 x 140 x 11.0",64.8,260.0,140.0,11.0,81.86,9.73,20.6,7017.0,2617.0,9.258,5.654,539.8,373.8,677.6,435.9,6169.0],  // MS derived; check availability
 ["260 x 140 x 12.5",72.3,260.0,140.0,12.5,92.07,8.2,17.8,7767.0,2876.0,9.184,5.589,597.4,410.9,755.7,484.5,6841.0],  // MS derived
 ["260 x 140 x 14.2",82.1,260.0,140.0,14.2,103.4,6.86,15.3,8560.0,3144.0,9.1,5.515,658.5,449.2,840.2,536.6,7555.0],  // MS derived; check availability
 ["260 x 140 x 16.0",90.3,260.0,140.0,16.0,115.0,5.75,13.2,9337.0,3400.0,9.01,5.437,718.3,485.8,925.1,588.3,8257.0],  // MS derived
 ["260 x 180 x 8.0",53.0,260.0,180.0,8.0,67.15,19.5,29.5,6390.0,3608.0,9.755,7.33,491.5,400.9,591.9,459.0,7221.0],  // MS derived; check availability
 ["260 x 180 x 8.8",58.1,260.0,180.0,8.8,73.51,17.5,26.5,6942.0,3912.0,9.717,7.295,534.0,434.7,645.4,500.0,7863.0],  // MS derived; check availability
 ["260 x 180 x 10.0",65.6,260.0,180.0,10.0,82.93,15.0,23.0,7741.0,4351.0,9.662,7.243,595.5,483.4,723.6,559.9,8798.0],  // MS derived; check availability
 ["260 x 180 x 11.0",71.7,260.0,180.0,11.0,90.66,13.4,20.6,8382.0,4700.0,9.615,7.2,644.7,522.2,787.1,608.4,9553.0],  // MS derived; check availability
 ["260 x 180 x 12.5",80.9,260.0,180.0,12.5,102.1,11.4,17.8,9299.0,5196.0,9.545,7.135,715.3,577.3,879.4,678.7,10640.0],  // MS derived; check availability
 ["260 x 180 x 14.2",91.0,260.0,180.0,14.2,114.7,9.68,15.3,10280.0,5719.0,9.465,7.06,790.6,635.5,979.8,754.7,11820.0],  // MS derived; check availability
 ["260 x 180 x 16.0",102.0,260.0,180.0,16.0,127.8,8.25,13.2,11250.0,6231.0,9.38,6.982,865.0,692.3,1081.0,831.1,12990.0],  // MS derived; check availability
 ["300 x 100 x 5.0",30.4,300.0,100.0,5.0,38.7,17.0,57.0,4150.0,731.0,10.3,4.34,276.0,146.0,354.0,161.0,2040.0],
 ["300 x 100 x 5.6",34.1,300.0,100.0,5.6,43.21,14.9,50.6,4598.0,805.8,10.32,4.318,306.5,161.2,393.3,178.9,2257.0],  // MS derived; check availability
 ["300 x 100 x 6.3",38.0,300.0,100.0,6.3,48.4,12.9,44.6,5110.0,890.0,10.3,4.29,341.0,178.0,439.0,199.0,2500.0],
 ["300 x 100 x 7.1",42.8,300.0,100.0,7.1,54.24,11.1,39.3,5683.0,980.9,10.24,4.253,378.9,196.2,489.7,221.2,2775.0],  // MS derived; check availability
 ["300 x 100 x 8.0",47.7,300.0,100.0,8.0,60.8,9.5,34.5,6300.0,1080.0,10.2,4.21,420.0,216.0,546.0,245.0,3070.0],
 ["300 x 100 x 8.8",52.5,300.0,100.0,8.8,66.47,8.36,31.1,6841.0,1160.0,10.14,4.178,456.0,232.1,594.4,266.2,3319.0],  // MS derived; check availability
 ["300 x 100 x 10.0",58.8,300.0,100.0,10.0,74.9,7.0,27.0,7610.0,1280.0,10.1,4.13,508.0,255.0,666.0,296.0,3680.0],
 ["300 x 100 x 11.0",64.8,300.0,100.0,11.0,81.86,6.09,24.3,8228.0,1364.0,10.03,4.082,548.5,272.8,723.0,320.1,3957.0],  // MS derived; check availability
 ["300 x 100 x 12.5",72.3,300.0,100.0,12.5,92.1,5.0,21.0,9100.0,1490.0,9.94,4.02,607.0,297.0,806.0,354.0,4350.0],
 ["300 x 100 x 14.2",81.1,300.0,100.0,14.2,103.0,4.04,18.1,10000.0,1610.0,9.85,3.94,669.0,321.0,896.0,390.0,4760.0],
 ["300 x 100 x 16.0",90.3,300.0,100.0,16.0,115.0,3.25,15.8,10900.0,1720.0,9.75,3.87,729.0,344.0,986.0,425.0,5140.0],
 ["300 x 100 x 17.5",97.7,300.0,100.0,17.5,124.0,2.71,14.1,11600.0,1800.0,9.66,3.8,775.0,360.0,1060.0,452.0,5420.0],
 ["300 x 200 x 5.0",38.3,300.0,200.0,5.0,48.7,37.0,57.0,6320.0,3400.0,11.4,8.35,421.0,340.0,501.0,380.0,6820.0],
 ["300 x 200 x 5.6",42.9,300.0,200.0,5.6,54.41,32.7,50.6,7025.0,3769.0,11.36,8.322,468.3,376.9,558.1,423.0,7593.0],  // MS derived; check availability
 ["300 x 200 x 6.0",45.7,300.0,200.0,6.0,58.17,30.3,47.0,7486.0,4013.0,11.34,8.305,499.1,401.3,595.8,451.3,8100.0],  // MS derived
 ["300 x 200 x 6.3",47.9,300.0,200.0,6.3,61.0,28.7,44.6,7830.0,4190.0,11.3,8.29,522.0,419.0,624.0,472.0,8480.0],
 ["300 x 200 x 7.1",54.0,300.0,200.0,7.1,68.44,25.2,39.3,8729.0,4667.0,11.29,8.258,581.9,466.7,697.6,527.9,9469.0],  // MS derived; check availability
 ["300 x 200 x 8.0",60.3,300.0,200.0,8.0,76.8,22.0,34.5,9720.0,5180.0,11.3,8.22,648.0,518.0,779.0,589.0,10600.0],
 ["300 x 200 x 8.8",66.4,300.0,200.0,8.8,84.07,19.7,31.1,10570.0,5631.0,11.21,8.184,704.9,563.1,850.7,642.6,11510.0],  // MS derived; check availability
 ["300 x 200 x 10.0",74.5,300.0,200.0,10.0,94.9,17.0,27.0,11800.0,6280.0,11.2,8.13,788.0,628.0,956.0,721.0,12900.0],
 ["300 x 200 x 11.0",82.1,300.0,200.0,11.0,103.9,15.2,24.3,12820.0,6795.0,11.11,8.089,854.9,679.5,1041.0,784.4,14040.0],  // MS derived; check availability
 ["300 x 200 x 12.0",88.5,300.0,200.0,12.0,112.7,13.7,22.0,13800.0,7294.0,11.06,8.045,919.8,729.4,1124.0,846.5,15140.0],  // MS derived
 ["300 x 200 x 12.5",91.9,300.0,200.0,12.5,117.0,13.0,21.0,14300.0,7540.0,11.0,8.02,952.0,754.0,1160.0,877.0,15700.0],
 ["300 x 200 x 14.2",103.0,300.0,200.0,14.2,132.0,11.1,18.1,15800.0,8330.0,11.0,7.95,1060.0,833.0,1300.0,978.0,17500.0],
 ["300 x 200 x 16.0",115.0,300.0,200.0,16.0,147.0,9.5,15.8,17400.0,9110.0,10.9,7.87,1160.0,911.0,1440.0,1080.0,19300.0],
 ["300 x 200 x 17.5",125.0,300.0,200.0,17.5,159.0,8.43,14.1,18600.0,9720.0,10.8,7.81,1240.0,972.0,1550.0,1160.0,20700.0],
 ["300 x 250 x 6.3",52.8,300.0,250.0,6.3,67.29,36.7,44.6,9188.0,6950.0,11.69,10.16,612.5,556.0,716.3,632.7,12160.0],  // MS derived; check availability
 ["300 x 250 x 7.1",59.5,300.0,250.0,7.1,75.54,32.2,39.3,10250.0,7749.0,11.65,10.13,683.5,619.9,801.6,707.9,13590.0],  // MS derived; check availability
 ["300 x 250 x 8.0",66.5,300.0,250.0,8.0,84.75,28.2,34.5,11420.0,8627.0,11.61,10.09,761.5,690.1,896.1,791.0,15190.0],  // MS derived; check availability
 ["300 x 250 x 8.8",73.3,300.0,250.0,8.8,92.87,25.4,31.1,12440.0,9388.0,11.57,10.05,829.3,751.0,978.8,863.7,16580.0],  // MS derived; check availability
 ["300 x 250 x 10.0",82.4,300.0,250.0,10.0,104.9,22.0,27.0,13920.0,10500.0,11.52,10.0,928.2,839.7,1101.0,970.7,18620.0],  // MS derived; check availability
 ["300 x 250 x 11.0",90.7,300.0,250.0,11.0,114.9,19.7,24.3,15120.0,11390.0,11.47,9.958,1008.0,911.2,1200.0,1058.0,20280.0],  // MS derived; check availability
 ["300 x 250 x 12.5",102.0,300.0,250.0,12.5,129.6,17.0,21.0,16860.0,12680.0,11.41,9.892,1124.0,1014.0,1345.0,1185.0,22710.0],  // MS derived; check availability
 ["300 x 250 x 14.2",116.0,300.0,250.0,14.2,146.0,14.6,18.1,18730.0,14070.0,11.33,9.818,1249.0,1126.0,1505.0,1325.0,25360.0],  // MS derived; check availability
 ["300 x 250 x 16.0",128.0,300.0,250.0,16.0,163.0,12.6,15.8,20620.0,15460.0,11.25,9.739,1375.0,1237.0,1668.0,1467.0,28060.0],  // MS derived; check availability
 ["350 x 150 x 5.0",38.3,350.0,150.0,5.0,48.7,27.0,67.0,7660.0,2050.0,12.5,6.49,437.0,274.0,543.0,301.0,5160.0],
 ["350 x 150 x 5.6",42.9,350.0,150.0,5.6,54.41,23.8,59.5,8506.0,2274.0,12.5,6.465,486.1,303.2,604.7,334.4,5733.0],  // MS derived; check availability
 ["350 x 150 x 6.3",47.9,350.0,150.0,6.3,61.0,20.8,52.6,9480.0,2520.0,12.5,6.43,542.0,337.0,676.0,373.0,6390.0],
 ["350 x 150 x 7.1",54.0,350.0,150.0,7.1,68.44,18.1,46.3,10570.0,2803.0,12.43,6.4,604.1,373.8,755.9,416.4,7122.0],  // MS derived; check availability
 ["350 x 150 x 8.0",60.3,350.0,150.0,8.0,76.8,15.8,40.8,11800.0,3100.0,12.4,6.36,673.0,414.0,844.0,464.0,7930.0],
 ["350 x 150 x 8.8",66.4,350.0,150.0,8.8,84.07,14.0,36.8,12810.0,3364.0,12.34,6.326,731.9,448.6,921.7,505.5,8622.0],  // MS derived; check availability
 ["350 x 150 x 10.0",74.5,350.0,150.0,10.0,94.9,12.0,32.0,14300.0,3740.0,12.3,6.27,818.0,498.0,1040.0,566.0,9630.0],
 ["350 x 150 x 11.0",82.1,350.0,150.0,11.0,103.9,10.6,28.8,15540.0,4032.0,12.23,6.231,887.9,537.6,1128.0,615.0,10450.0],  // MS derived; check availability
 ["350 x 150 x 12.5",91.9,350.0,150.0,12.5,117.0,9.0,25.0,17300.0,4450.0,12.2,6.17,988.0,593.0,1260.0,686.0,11600.0],
 ["350 x 150 x 14.2",103.0,350.0,150.0,14.2,132.0,7.56,21.6,19200.0,4890.0,12.1,6.09,1100.0,652.0,1410.0,763.0,12900.0],
 ["350 x 150 x 16.0",115.0,350.0,150.0,16.0,147.0,6.38,18.9,21100.0,5320.0,12.0,6.01,1200.0,709.0,1560.0,840.0,14100.0],
 ["350 x 150 x 17.5",125.0,350.0,150.0,17.5,159.0,5.57,17.0,22600.0,5640.0,11.9,5.95,1290.0,752.0,1680.0,900.0,15100.0],
 ["350 x 200 x 6.3",53.0,350.0,200.0,6.3,67.29,28.7,52.6,11340.0,4785.0,12.98,8.433,648.1,478.5,784.1,533.4,10520.0],  // MS derived; check availability
 ["350 x 200 x 7.1",59.5,350.0,200.0,7.1,75.54,25.2,46.3,12660.0,5328.0,12.95,8.398,723.4,532.8,877.6,596.4,11750.0],  // MS derived; check availability
 ["350 x 200 x 8.0",66.8,350.0,200.0,8.0,84.75,22.0,40.8,14110.0,5922.0,12.9,8.359,806.3,592.2,981.2,666.0,13110.0],  // MS derived; check availability
 ["350 x 200 x 8.8",73.3,350.0,200.0,8.8,92.87,19.7,36.8,15370.0,6436.0,12.86,8.325,878.3,643.6,1072.0,726.7,14300.0],  // MS derived; check availability
 ["350 x 200 x 10.0",82.8,350.0,200.0,10.0,104.9,17.0,32.0,17210.0,7181.0,12.81,8.273,983.5,718.1,1205.0,815.9,16040.0],  // MS derived; check availability
 ["350 x 200 x 11.0",90.7,350.0,200.0,11.0,114.9,15.2,28.8,18700.0,7779.0,12.76,8.229,1069.0,777.9,1314.0,888.4,17450.0],  // MS derived; check availability
 ["350 x 200 x 12.5",102.0,350.0,200.0,12.5,129.6,13.0,25.0,20860.0,8637.0,12.69,8.165,1192.0,863.7,1474.0,994.2,19510.0],  // MS derived; check availability
 ["350 x 200 x 14.2",116.0,350.0,200.0,14.2,146.0,11.1,21.6,23190.0,9556.0,12.61,8.091,1325.0,955.6,1649.0,1110.0,21740.0],  // MS derived; check availability
 ["350 x 200 x 16.0",129.0,350.0,200.0,16.0,163.0,9.5,18.9,25550.0,10470.0,12.52,8.013,1460.0,1047.0,1828.0,1227.0,24000.0],  // MS derived; check availability
 ["350 x 250 x 6.3",57.8,350.0,250.0,6.3,73.6,36.7,52.6,13200.0,7880.0,13.4,10.4,754.0,631.0,892.0,709.0,15200.0],
 ["350 x 250 x 7.1",65.1,350.0,250.0,7.1,82.64,32.2,46.3,14750.0,8796.0,13.36,10.32,842.7,703.7,999.3,794.1,17020.0],  // MS derived; check availability
 ["350 x 250 x 8.0",72.8,350.0,250.0,8.0,92.8,28.3,40.8,16400.0,9800.0,13.3,10.3,940.0,784.0,1120.0,888.0,19000.0],
 ["350 x 250 x 8.8",80.2,350.0,250.0,8.8,101.7,25.4,36.8,17930.0,10670.0,13.28,10.24,1025.0,853.4,1222.0,969.9,20780.0],  // MS derived; check availability
 ["350 x 250 x 10.0",90.2,350.0,250.0,10.0,115.0,22.0,32.0,20100.0,11900.0,13.2,10.2,1150.0,955.0,1380.0,1090.0,23400.0],
 ["350 x 250 x 11.0",99.4,350.0,250.0,11.0,125.9,19.7,28.8,21860.0,12960.0,13.18,10.15,1249.0,1037.0,1501.0,1189.0,25450.0],  // MS derived; check availability
 ["350 x 250 x 12.5",112.0,350.0,250.0,12.5,142.0,17.0,25.0,24400.0,14400.0,13.1,10.1,1400.0,1160.0,1680.0,1330.0,28500.0],
 ["350 x 250 x 14.2",126.0,350.0,250.0,14.2,160.0,14.6,21.6,27200.0,16000.0,13.0,10.0,1550.0,1280.0,1890.0,1490.0,31900.0],
 ["350 x 250 x 16.0",141.0,350.0,250.0,16.0,179.0,12.6,18.9,30000.0,17700.0,12.9,9.93,1720.0,1410.0,2100.0,1660.0,35300.0],
 ["350 x 250 x 17.5",153.0,350.0,250.0,17.5,194.0,11.3,17.0,32200.0,18900.0,12.9,9.87,1840.0,1510.0,2260.0,1780.0,38100.0],
 ["400 x 150 x 6.3",52.8,400.0,150.0,6.3,67.29,20.8,60.5,13250.0,2851.0,14.04,6.509,662.7,380.1,836.2,418.3,7595.0],  // MS derived; check availability
 ["400 x 150 x 7.1",59.5,400.0,150.0,7.1,75.54,18.1,53.3,14790.0,3166.0,13.99,6.474,739.7,422.1,935.9,467.1,8467.0],  // MS derived; check availability
 ["400 x 150 x 8.0",66.5,400.0,150.0,8.0,84.75,15.8,47.0,16490.0,3509.0,13.95,6.435,824.4,467.9,1046.0,520.9,9424.0],  // MS derived; check availability
 ["400 x 150 x 8.8",73.3,400.0,150.0,8.8,92.87,14.0,42.5,17960.0,3804.0,13.91,6.4,898.0,507.2,1143.0,567.6,10250.0],  // MS derived; check availability
 ["400 x 150 x 10.0",82.4,400.0,150.0,10.0,104.9,12.0,37.0,20110.0,4227.0,13.84,6.347,1006.0,563.7,1285.0,636.0,11460.0],  // MS derived; check availability
 ["400 x 150 x 11.0",90.7,400.0,150.0,11.0,114.9,10.6,33.4,21850.0,4564.0,13.79,6.304,1092.0,608.6,1401.0,691.4,12430.0],  // MS derived; check availability
 ["400 x 150 x 12.5",102.0,400.0,150.0,12.5,129.6,9.0,29.0,24370.0,5043.0,13.71,6.238,1218.0,672.4,1571.0,771.8,13830.0],  // MS derived; check availability
 ["400 x 150 x 14.2",116.0,400.0,150.0,14.2,146.0,7.56,25.2,27100.0,5547.0,13.62,6.164,1355.0,739.6,1758.0,858.9,15330.0],  // MS derived; check availability
 ["400 x 150 x 16.0",128.0,400.0,150.0,16.0,163.0,6.38,22.0,29840.0,6038.0,13.53,6.086,1492.0,805.1,1948.0,946.7,16830.0],  // MS derived; check availability
 ["400 x 200 x 6.0",55.1,400.0,200.0,6.0,70.17,30.3,63.7,15000.0,5142.0,14.62,8.56,749.9,514.2,916.6,567.7,12050.0],  // MS derived
 ["400 x 200 x 6.3",57.8,400.0,200.0,6.3,73.6,28.7,60.5,15700.0,5380.0,14.6,8.55,785.0,538.0,960.0,594.0,12600.0],
 ["400 x 200 x 7.1",65.1,400.0,200.0,7.1,82.64,25.2,53.3,17530.0,5989.0,14.57,8.513,876.7,598.9,1075.0,664.8,14100.0],  // MS derived; check availability
 ["400 x 200 x 8.0",72.8,400.0,200.0,8.0,92.8,22.0,47.0,19600.0,6660.0,14.5,8.47,978.0,666.0,1200.0,743.0,15700.0],
 ["400 x 200 x 8.8",80.2,400.0,200.0,8.8,101.7,19.7,42.5,21330.0,7241.0,14.48,8.439,1066.0,724.1,1315.0,810.8,17160.0],  // MS derived; check availability
 ["400 x 200 x 10.0",90.2,400.0,200.0,10.0,115.0,17.0,37.0,23900.0,8080.0,14.4,8.39,1200.0,808.0,1480.0,911.0,19300.0],
 ["400 x 200 x 11.0",99.4,400.0,200.0,11.0,125.9,15.2,33.4,26010.0,8762.0,14.38,8.344,1301.0,876.2,1615.0,992.3,20960.0],  // MS derived; check availability
 ["400 x 200 x 12.0",107.0,400.0,200.0,12.0,136.7,13.7,30.3,28060.0,9418.0,14.33,8.3,1403.0,941.8,1748.0,1072.0,22620.0],  // MS derived
 ["400 x 200 x 12.5",112.0,400.0,200.0,12.5,142.0,13.0,29.0,29100.0,9740.0,14.3,8.28,1450.0,974.0,1810.0,1110.0,23400.0],
 ["400 x 200 x 14.2",126.0,400.0,200.0,14.2,160.0,11.1,25.2,32400.0,10800.0,14.2,8.21,1620.0,1080.0,2030.0,1240.0,26100.0],
 ["400 x 200 x 16.0",141.0,400.0,200.0,16.0,179.0,9.5,22.0,35700.0,11800.0,14.1,8.13,1790.0,1180.0,2260.0,1370.0,28900.0],
 ["400 x 200 x 17.5",153.0,400.0,200.0,17.5,194.0,8.43,19.9,38400.0,12600.0,14.1,8.06,1920.0,1260.0,2440.0,1480.0,31100.0],
 ["400 x 300 x 8.0",85.4,400.0,300.0,8.0,109.0,34.5,47.0,25700.0,16500.0,15.4,12.3,1280.0,1100.0,1520.0,1250.0,31000.0],
 ["400 x 300 x 8.8",94.0,400.0,300.0,8.8,119.3,31.1,42.5,28060.0,18040.0,15.34,12.3,1403.0,1202.0,1659.0,1363.0,33910.0],  // MS derived; check availability
 ["400 x 300 x 10.0",106.0,400.0,300.0,10.0,135.0,27.0,37.0,31500.0,20200.0,15.3,12.2,1580.0,1350.0,1870.0,1540.0,38200.0],
 ["400 x 300 x 11.0",117.0,400.0,300.0,11.0,147.9,24.3,33.4,34340.0,22020.0,15.24,12.2,1717.0,1468.0,2043.0,1677.0,41680.0],  // MS derived; check availability
 ["400 x 300 x 12.5",131.0,400.0,300.0,12.5,167.0,21.0,29.0,38500.0,24600.0,15.2,12.1,1920.0,1640.0,2300.0,1880.0,46800.0],
 ["400 x 300 x 14.2",148.0,400.0,300.0,14.2,189.0,18.1,25.2,43000.0,27400.0,15.1,12.1,2150.0,1830.0,2580.0,2110.0,52500.0],
 ["400 x 300 x 16.0",166.0,400.0,300.0,16.0,211.0,15.8,22.0,47500.0,30300.0,15.0,12.0,2380.0,2020.0,2870.0,2350.0,58300.0],
 ["400 x 300 x 17.5",180.0,400.0,300.0,17.5,229.0,14.1,19.9,51200.0,32600.0,14.9,11.9,2560.0,2170.0,3110.0,2540.0,63000.0],
 ["450 x 150 x 8.0",73.1,450.0,150.0,8.0,92.75,15.8,53.2,22270.0,3913.0,15.49,6.495,989.6,521.7,1268.0,577.7,10940.0],  // MS derived; check availability
 ["450 x 150 x 8.8",80.2,450.0,150.0,8.8,101.7,14.0,48.1,24270.0,4243.0,15.45,6.46,1079.0,565.7,1386.0,629.8,11910.0],  // MS derived; check availability
 ["450 x 150 x 10.0",90.7,450.0,150.0,10.0,114.9,12.0,42.0,27210.0,4718.0,15.39,6.407,1209.0,629.1,1560.0,706.0,13310.0],  // MS derived; check availability
 ["450 x 150 x 11.0",99.4,450.0,150.0,11.0,125.9,10.6,37.9,29600.0,5097.0,15.33,6.364,1315.0,679.6,1702.0,767.9,14440.0],  // MS derived; check availability
 ["450 x 150 x 12.5",112.0,450.0,150.0,12.5,142.1,9.0,33.0,33060.0,5635.0,15.25,6.298,1469.0,751.4,1911.0,857.7,16070.0],  // MS derived; check availability
 ["450 x 150 x 14.2",127.0,450.0,150.0,14.2,160.2,7.56,28.7,36830.0,6204.0,15.16,6.224,1637.0,827.2,2140.0,955.4,17820.0],  // MS derived; check availability
 ["450 x 150 x 16.0",142.0,450.0,150.0,16.0,179.0,6.38,25.1,40630.0,6760.0,15.07,6.145,1806.0,901.3,2376.0,1054.0,19570.0],  // MS derived; check availability
 ["450 x 250 x 8.0",85.4,450.0,250.0,8.0,109.0,28.3,53.3,30100.0,12100.0,16.6,10.6,1340.0,971.0,1620.0,1080.0,27100.0],
 ["450 x 250 x 8.8",94.0,450.0,250.0,8.8,119.3,25.4,48.1,32840.0,13230.0,16.59,10.53,1460.0,1058.0,1774.0,1182.0,29590.0],  // MS derived; check availability
 ["450 x 250 x 10.0",106.0,450.0,250.0,10.0,135.0,22.0,42.0,36900.0,14800.0,16.5,10.5,1640.0,1180.0,2000.0,1330.0,33300.0],
 ["450 x 250 x 11.0",117.0,450.0,250.0,11.0,147.9,19.7,37.9,40200.0,16110.0,16.49,10.44,1787.0,1288.0,2185.0,1452.0,36300.0],  // MS derived; check availability
 ["450 x 250 x 12.0",126.0,450.0,250.0,12.0,160.7,17.8,34.5,43430.0,17360.0,16.44,10.39,1930.0,1389.0,2367.0,1572.0,39260.0],  // MS derived
 ["450 x 250 x 12.5",131.0,450.0,250.0,12.5,167.0,17.0,33.0,45000.0,18000.0,16.4,10.4,2000.0,1440.0,2460.0,1630.0,40700.0],
 ["450 x 250 x 14.2",148.0,450.0,250.0,14.2,189.0,14.6,28.7,50300.0,20000.0,16.3,10.3,2240.0,1600.0,2760.0,1830.0,45600.0],
 ["450 x 250 x 16.0",166.0,450.0,250.0,16.0,211.0,12.6,25.1,55700.0,22000.0,16.2,10.2,2480.0,1760.0,3070.0,2030.0,50500.0],
 ["450 x 250 x 17.5",180.0,450.0,250.0,17.5,229.0,11.3,22.7,60000.0,23700.0,16.2,10.2,2670.0,1890.0,3320.0,2190.0,54600.0],
 ["500 x 200 x 8.0",85.4,500.0,200.0,8.0,109.0,22.0,59.5,34000.0,8140.0,17.7,8.65,1360.0,814.0,1710.0,896.0,21100.0],
 ["500 x 200 x 8.8",94.0,500.0,200.0,8.8,119.3,19.7,53.8,37170.0,8850.0,17.65,8.614,1487.0,885.0,1867.0,979.1,23050.0],  // MS derived; check availability
 ["500 x 200 x 10.0",106.0,500.0,200.0,10.0,135.0,17.0,47.0,41800.0,9890.0,17.6,8.56,1670.0,989.0,2100.0,1100.0,25900.0],
 ["500 x 200 x 11.0",117.0,500.0,200.0,11.0,147.9,15.2,42.5,45490.0,10730.0,17.54,8.518,1820.0,1073.0,2299.0,1200.0,28170.0],  // MS derived; check availability
 ["500 x 200 x 12.5",131.0,500.0,200.0,12.5,167.0,13.0,37.0,51000.0,11900.0,17.5,8.45,2040.0,1190.0,2590.0,1350.0,31500.0],
 ["500 x 200 x 14.2",148.0,500.0,200.0,14.2,189.0,11.1,32.2,56900.0,13200.0,17.4,8.38,2280.0,1320.0,2900.0,1500.0,35200.0],
 ["500 x 200 x 16.0",166.0,500.0,200.0,16.0,211.0,9.5,28.3,63000.0,14500.0,17.3,8.3,2520.0,1450.0,3230.0,1670.0,38900.0],
 ["500 x 200 x 17.5",180.0,500.0,200.0,17.5,229.0,8.43,25.6,67900.0,15600.0,17.2,8.24,2720.0,1560.0,3500.0,1800.0,41800.0],
 ["500 x 300 x 8.0",97.9,500.0,300.0,8.0,125.0,34.5,59.5,43700.0,20000.0,18.7,12.6,1750.0,1330.0,2100.0,1480.0,42600.0],
 ["500 x 300 x 8.8",108.0,500.0,300.0,8.8,136.9,31.1,53.8,47780.0,21770.0,18.68,12.61,1911.0,1451.0,2300.0,1619.0,46550.0],  // MS derived; check availability
 ["500 x 300 x 10.0",122.0,500.0,300.0,10.0,155.0,27.0,47.0,53800.0,24400.0,18.6,12.6,2150.0,1630.0,2600.0,1830.0,52400.0],
 ["500 x 300 x 11.0",134.0,500.0,300.0,11.0,169.9,24.3,42.5,58650.0,26610.0,18.58,12.52,2346.0,1774.0,2837.0,1995.0,57280.0],  // MS derived; check availability
 ["500 x 300 x 12.0",145.0,500.0,300.0,12.0,184.7,22.0,38.7,63450.0,28740.0,18.53,12.47,2538.0,1916.0,3077.0,2161.0,62040.0],  // MS derived
 ["500 x 300 x 12.5",151.0,500.0,300.0,12.5,192.0,21.0,37.0,65800.0,29800.0,18.5,12.5,2630.0,1980.0,3200.0,2240.0,64400.0],
 ["500 x 300 x 14.2",170.0,500.0,300.0,14.2,217.0,18.1,32.2,73700.0,33200.0,18.4,12.4,2950.0,2220.0,3590.0,2520.0,72200.0],
 ["500 x 300 x 16.0",191.0,500.0,300.0,16.0,243.0,15.8,28.3,81800.0,36800.0,18.3,12.3,3270.0,2450.0,4000.0,2800.0,80300.0],
 ["500 x 300 x 17.5",208.0,500.0,300.0,17.5,264.0,14.1,25.6,88300.0,39600.0,18.3,12.2,3530.0,2640.0,4340.0,3040.0,86900.0],
 ["500 x 300 x 20.0",235.0,500.0,300.0,20.0,299.7,12.0,22.0,98780.0,44080.0,18.15,12.13,3951.0,2939.0,4885.0,3408.0,97450.0]  // MS derived
].map(a=>({key:a[0],mass:a[1],D:a[2],B:a[3],t:a[4],A:a[5],bT:a[6],dt:a[7],
  Ix:a[8],Iy:a[9],rx:a[10],ry:a[11],Zx:a[12],Zy:a[13],Sx:a[14],Sy:a[15],J:a[16]}));

const RHSmap=Object.fromEntries(RHS.map(s=>[s.key,s]));

/* Cold-formed RHS (EN 10219, S355) from the reference software 2025 UK library. Same row layout as RHS above.
   DATA ONLY: activeSectionBase() in js/02-section-data.js has no cold-formed RHS switch yet (the RHS family is
   normalised as boxType "HF"), so these rows are not reachable from the UI until an rhsType selector mirroring
   shsType is wired in js/02-section-data.js, js/03-state-ui.js, js/07-wiring.js and index.html. When wired, use
   boxType "CF" (BS 5950 curve c, EN 10219 S355 design strength, EC3 c = h - 3t per P363 as SHS_CF does).
   bT = (B - 5t)/t, dt = (D - 5t)/t (BS 5950 Table 11 cold-formed flat widths, matching SHS_CF).
   A, Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, J derived with EN 10219-2:2006 Annex B (ro = 2.0t / 2.5t / 3.0t for
   t <= 6 / 6 < t <= 10 / t > 10, ri = ro - t; J per the ctBoxEN10210() convention with mean corner radius (ro + ri)/2). */
const RHS_CF=[
 // key, mass, D(h), B(b), t, A, bT(flange cf/t), dt(web cw/t), Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, J
 ["50 x 25 x 2.0",2.15,50.0,25.0,2.0,2.737,7.5,20.0,8.384,2.809,1.75,1.013,3.353,2.247,4.262,2.616,7.062],  // MS derived; check availability
 ["50 x 25 x 2.5",2.62,50.0,25.0,2.5,3.339,5.0,15.0,9.887,3.278,1.721,0.9909,3.955,2.623,5.105,3.119,8.43],  // MS derived
 ["50 x 25 x 3.0",3.07,50.0,25.0,3.0,3.908,3.33,11.7,11.17,3.667,1.691,0.9687,4.469,2.934,5.863,3.565,9.637],  // MS derived
 ["50 x 30 x 2.0",2.31,50.0,30.0,2.0,2.937,10.0,20.0,9.536,4.293,1.802,1.209,3.815,2.862,4.742,3.325,9.767],  // MS derived; check availability
 ["50 x 30 x 2.5",2.82,50.0,30.0,2.5,3.589,7.0,15.0,11.3,5.052,1.774,1.186,4.519,3.368,5.699,3.985,11.74],  // MS derived
 ["50 x 30 x 3.0",3.3,50.0,30.0,3.0,4.208,5.0,11.7,12.83,5.7,1.746,1.164,5.132,3.8,6.568,4.579,13.53],  // MS derived
 ["50 x 30 x 4.0",4.2,50.0,30.0,4.0,5.348,2.5,7.5,15.25,6.693,1.689,1.119,6.1,4.462,8.049,5.581,16.53],  // MS derived
 ["60 x 30 x 3.0",3.77,60.0,30.0,3.0,4.808,5.0,15.0,20.5,6.798,2.065,1.189,6.834,4.532,8.822,5.389,17.48],  // MS derived
 ["60 x 30 x 4.0",4.83,60.0,30.0,4.0,6.148,2.5,10.0,24.7,8.055,2.005,1.145,8.234,5.37,10.92,6.621,21.47],  // MS derived
 ["60 x 40 x 2.5",3.6,60.0,40.0,2.5,4.589,11.0,19.0,22.07,11.74,2.193,1.599,7.357,5.868,9.056,6.842,25.14],  // MS derived; check availability
 ["60 x 40 x 3.0",4.25,60.0,40.0,3.0,5.408,8.33,15.0,25.38,13.44,2.166,1.576,8.46,6.72,10.53,7.944,29.28],  // MS derived
 ["60 x 40 x 4.0",5.45,60.0,40.0,4.0,6.948,5.0,10.0,30.99,16.28,2.112,1.531,10.33,8.14,13.16,9.895,36.67],  // MS derived
 ["60 x 40 x 5.0",6.56,60.0,40.0,5.0,8.356,3.0,7.0,35.33,18.43,2.056,1.485,11.78,9.213,15.38,11.52,42.85],  // MS derived
 ["70 x 40 x 3.0",4.72,70.0,40.0,3.0,6.008,8.33,18.3,37.31,15.5,2.492,1.606,10.66,7.749,13.39,9.054,36.49],  // MS derived; check availability
 ["70 x 40 x 4.0",6.08,70.0,40.0,4.0,7.748,5.0,12.5,45.95,18.88,2.435,1.561,13.13,9.441,16.84,11.33,45.84],  // MS derived; check availability
 ["70 x 40 x 5.0",7.34,70.0,40.0,5.0,9.356,3.0,9.0,52.88,21.51,2.377,1.516,15.11,10.75,19.81,13.27,53.77],  // MS derived; check availability
 ["70 x 50 x 3.0",5.19,70.0,50.0,3.0,6.608,11.7,18.3,44.05,26.1,2.582,1.987,12.59,10.44,15.4,12.21,53.62],  // MS derived; check availability
 ["70 x 50 x 4.0",6.71,70.0,50.0,4.0,8.548,7.5,12.5,54.67,32.22,2.529,1.942,15.62,12.89,19.48,15.41,68.07],  // MS derived
 ["70 x 50 x 5.0",8.13,70.0,50.0,5.0,10.36,5.0,9.0,63.46,37.2,2.475,1.895,18.13,14.88,23.06,18.2,80.77],  // MS derived
 ["80 x 40 x 3.0",5.19,80.0,40.0,3.0,6.608,8.33,21.7,52.25,17.56,2.812,1.63,13.06,8.778,16.54,10.16,43.88],  // MS derived
 ["80 x 40 x 4.0",6.71,80.0,40.0,4.0,8.548,5.0,15.0,64.79,21.49,2.753,1.585,16.2,10.74,20.91,12.77,55.24],  // MS derived
 ["80 x 40 x 5.0",8.13,80.0,40.0,5.0,10.36,3.0,11.0,75.11,24.59,2.693,1.541,18.78,12.3,24.74,15.02,64.97],  // MS derived
 ["80 x 50 x 3.0",5.66,80.0,50.0,3.0,7.208,11.7,21.7,61.15,29.42,2.913,2.02,15.29,11.77,18.85,13.62,65.0],  // MS derived; check availability
 ["80 x 50 x 4.0",7.34,80.0,50.0,4.0,9.348,7.5,15.0,76.36,36.46,2.858,1.975,19.09,14.59,23.95,17.25,82.7],  // MS derived
 ["80 x 50 x 5.0",8.91,80.0,50.0,5.0,11.36,5.0,11.0,89.19,42.29,2.803,1.93,22.3,16.92,28.49,20.45,98.4],  // MS derived
 ["80 x 60 x 3.0",6.13,80.0,60.0,3.0,7.808,15.0,21.7,70.05,44.89,2.995,2.398,17.51,14.96,21.16,17.37,88.35],  // MS derived
 ["80 x 60 x 3.5",7.06,80.0,60.0,3.5,8.995,12.1,17.9,79.3,50.72,2.969,2.375,19.83,16.91,24.15,19.81,101.0],  // MS derived; check availability
 ["80 x 60 x 4.0",7.97,80.0,60.0,4.0,10.15,10.0,15.0,87.92,56.12,2.943,2.352,21.98,18.71,26.99,22.12,113.1],  // MS derived
 ["80 x 60 x 5.0",9.7,80.0,60.0,5.0,12.36,7.0,11.0,103.3,65.66,2.891,2.305,25.82,21.89,32.24,26.38,135.5],  // MS derived
 ["90 x 50 x 3.0",6.13,90.0,50.0,3.0,7.808,11.7,25.0,81.85,32.74,3.238,2.048,18.19,13.1,22.6,15.03,76.67],  // MS derived; check availability
 ["90 x 50 x 3.5",7.06,90.0,50.0,3.5,8.995,9.29,20.7,92.66,36.89,3.21,2.025,20.59,14.76,25.79,17.12,87.47],  // MS derived
 ["90 x 50 x 4.0",7.97,90.0,50.0,4.0,10.15,7.5,17.5,102.7,40.71,3.181,2.003,22.82,16.28,28.82,19.09,97.7],  // MS derived; check availability
 ["90 x 50 x 5.0",9.7,90.0,50.0,5.0,12.36,5.0,13.0,120.6,47.37,3.124,1.958,26.8,18.95,34.41,22.7,116.5],  // MS derived; check availability
 ["100 x 40 x 3.0",6.13,100.0,40.0,3.0,7.808,8.33,28.3,92.34,21.67,3.439,1.666,18.47,10.84,23.75,12.38,59.05],  // MS derived; check availability
 ["100 x 40 x 4.0",7.97,100.0,40.0,4.0,10.15,5.0,20.0,115.7,26.69,3.377,1.622,23.14,13.35,30.26,15.65,74.53],  // MS derived; check availability
 ["100 x 40 x 5.0",9.7,100.0,40.0,5.0,12.36,3.0,15.0,135.6,30.76,3.313,1.578,27.12,15.38,36.09,18.52,87.92],  // MS derived; check availability
 ["100 x 50 x 3.0",6.6,100.0,50.0,3.0,8.408,11.7,28.3,106.5,36.06,3.558,2.071,21.29,14.42,26.66,16.44,88.56],  // MS derived
 ["100 x 50 x 4.0",8.59,100.0,50.0,4.0,10.95,7.5,20.0,134.1,44.95,3.5,2.026,26.83,17.98,34.1,20.93,113.0],  // MS derived
 ["100 x 50 x 5.0",10.5,100.0,50.0,5.0,13.36,5.0,15.0,158.2,52.45,3.441,1.982,31.64,20.98,40.84,24.95,134.9],  // MS derived
 ["100 x 50 x 6.0",12.3,100.0,50.0,6.0,15.63,3.33,11.7,178.8,58.67,3.381,1.937,35.75,23.47,46.9,28.52,154.2],  // MS derived
 ["100 x 60 x 3.0",7.07,100.0,60.0,3.0,9.008,15.0,28.3,120.6,54.65,3.659,2.463,24.11,18.22,29.57,20.79,121.7],  // MS derived
 ["100 x 60 x 3.5",8.16,100.0,60.0,3.5,10.39,12.1,23.6,137.1,61.91,3.631,2.44,27.41,20.64,33.84,23.77,139.3],  // MS derived; check availability
 ["100 x 60 x 4.0",9.22,100.0,60.0,4.0,11.75,10.0,20.0,152.6,68.68,3.604,2.418,30.52,22.89,37.94,26.6,156.3],  // MS derived
 ["100 x 60 x 5.0",11.3,100.0,60.0,5.0,14.36,7.0,15.0,180.8,80.83,3.548,2.373,36.15,26.94,45.59,31.88,187.9],  // MS derived
 ["100 x 60 x 6.0",13.2,100.0,60.0,6.0,16.83,5.0,11.7,205.3,91.2,3.492,2.328,41.06,30.4,52.54,36.64,216.4],  // MS derived
 ["100 x 80 x 3.0",8.01,100.0,80.0,3.0,10.21,21.7,28.3,148.8,105.6,3.818,3.217,29.76,26.41,35.39,30.4,196.1],  // MS derived; check availability
 ["100 x 80 x 4.0",10.5,100.0,80.0,4.0,13.35,15.0,20.0,189.5,134.2,3.768,3.17,37.89,33.54,45.62,39.15,253.8],  // MS derived
 ["100 x 80 x 5.0",12.8,100.0,80.0,5.0,16.36,11.0,15.0,225.9,159.6,3.717,3.124,45.19,39.9,55.09,47.24,307.5],  // MS derived
 ["100 x 80 x 6.0",14.9,100.0,80.0,6.0,19.23,8.33,11.7,258.4,182.1,3.665,3.077,51.68,45.53,63.82,54.67,357.4],  // MS derived
 ["120 x 40 x 3.0",7.07,120.0,40.0,3.0,9.008,8.33,35.0,148.0,25.79,4.054,1.692,24.67,12.89,32.16,14.6,74.56],  // MS derived; check availability
 ["120 x 40 x 4.0",9.22,120.0,40.0,4.0,11.75,5.0,25.0,186.9,31.9,3.989,1.648,31.15,15.95,41.21,18.53,94.23],  // MS derived; check availability
 ["120 x 40 x 5.0",11.3,120.0,40.0,5.0,14.36,3.0,19.0,220.8,36.93,3.922,1.604,36.8,18.46,49.45,22.02,111.4],  // MS derived; check availability
 ["120 x 60 x 3.0",8.01,120.0,60.0,3.0,10.21,15.0,35.0,189.1,64.4,4.304,2.512,31.52,21.47,39.18,24.21,156.3],  // MS derived; check availability
 ["120 x 60 x 3.5",9.26,120.0,60.0,3.5,11.79,12.1,29.3,215.6,73.09,4.276,2.489,35.94,24.36,44.94,27.72,179.2],  // MS derived
 ["120 x 60 x 4.0",10.5,120.0,60.0,4.0,13.35,10.0,25.0,240.7,81.25,4.247,2.467,40.12,27.08,50.49,31.08,201.1],  // MS derived
 ["120 x 60 x 5.0",12.8,120.0,60.0,5.0,16.36,7.0,19.0,287.0,95.99,4.189,2.423,47.83,32.0,60.95,37.38,242.2],  // MS derived
 ["120 x 60 x 6.0",15.1,120.0,60.0,6.0,19.23,5.0,15.0,328.0,108.8,4.13,2.378,54.67,36.26,70.57,43.12,279.7],  // MS derived
 ["120 x 80 x 3.0",8.96,120.0,80.0,3.0,11.41,21.7,35.0,230.2,123.4,4.492,3.289,38.37,30.86,46.2,35.02,255.5],  // MS derived; check availability
 ["120 x 80 x 4.0",11.7,120.0,80.0,4.0,14.95,15.0,25.0,294.6,157.3,4.439,3.244,49.1,39.32,59.77,45.23,331.2],  // MS derived
 ["120 x 80 x 5.0",14.4,120.0,80.0,5.0,18.36,11.0,19.0,353.1,187.8,4.386,3.198,58.86,46.94,72.45,54.74,402.3],  // MS derived
 ["120 x 80 x 6.0",17.0,120.0,80.0,6.0,21.63,8.33,15.0,406.1,215.0,4.332,3.153,67.68,53.76,84.25,63.55,468.5],  // MS derived
 ["120 x 80 x 8.0",21.4,120.0,80.0,8.0,27.24,5.0,10.0,475.8,251.7,4.179,3.039,79.31,62.92,102.0,76.93,584.0],  // MS derived
 ["140 x 80 x 3.0",9.9,140.0,80.0,3.0,12.61,21.7,41.7,334.4,141.2,5.15,3.347,47.77,35.31,58.2,39.64,317.1],  // MS derived; check availability
 ["140 x 80 x 4.0",13.0,140.0,80.0,4.0,16.55,15.0,30.0,429.6,180.4,5.095,3.302,61.37,45.1,75.51,51.31,411.6],  // MS derived; check availability
 ["140 x 80 x 5.0",16.0,140.0,80.0,5.0,20.36,11.0,23.0,517.1,215.9,5.04,3.257,73.87,53.99,91.8,62.24,500.5],  // MS derived
 ["140 x 80 x 6.0",18.9,140.0,80.0,6.0,24.03,8.33,18.3,597.0,248.0,4.984,3.212,85.29,61.99,107.1,72.43,583.8],  // MS derived
 ["140 x 80 x 8.0",23.9,140.0,80.0,8.0,30.44,5.0,12.5,708.1,293.3,4.823,3.104,101.2,73.33,130.8,88.45,731.4],  // MS derived
 ["140 x 80 x 10.0",28.7,140.0,80.0,10.0,36.57,3.0,9.0,803.7,330.5,4.688,3.006,114.8,82.62,152.4,102.7,851.0],  // MS derived
 ["150 x 100 x 3.0",11.3,150.0,100.0,3.0,14.41,28.3,45.0,460.6,247.6,5.654,4.146,61.42,49.53,73.48,55.76,507.2],  // MS derived; check availability
 ["150 x 100 x 4.0",14.9,150.0,100.0,4.0,18.95,20.0,32.5,594.6,318.6,5.602,4.1,79.28,63.71,95.67,72.5,661.6],  // MS derived
 ["150 x 100 x 5.0",18.3,150.0,100.0,5.0,23.36,15.0,25.0,719.2,384.0,5.549,4.055,95.89,76.8,116.7,88.34,808.7],  // MS derived
 ["150 x 100 x 6.0",21.7,150.0,100.0,6.0,27.63,11.7,20.0,834.7,444.2,5.496,4.009,111.3,88.84,136.7,103.3,948.3],  // MS derived
 ["150 x 100 x 8.0",27.7,150.0,100.0,8.0,35.24,7.5,13.8,1008.0,535.7,5.348,3.899,134.4,107.1,169.2,127.9,1206.0],  // MS derived
 ["150 x 100 x 10.0",33.4,150.0,100.0,10.0,42.57,5.0,10.0,1162.0,614.4,5.224,3.799,154.9,122.9,199.2,150.2,1426.0],  // MS derived
 ["160 x 80 x 3.0",10.8,160.0,80.0,3.0,13.81,21.7,48.3,463.8,159.0,5.796,3.394,57.98,39.76,71.41,44.26,380.3],  // MS derived; check availability
 ["160 x 80 x 4.0",14.2,160.0,80.0,4.0,18.15,15.0,35.0,597.7,203.5,5.739,3.349,74.71,50.89,92.86,57.39,494.1],  // MS derived
 ["160 x 80 x 5.0",17.5,160.0,80.0,5.0,22.36,11.0,27.0,721.7,244.1,5.682,3.304,90.21,61.03,113.2,69.74,601.3],  // MS derived
 ["160 x 80 x 6.0",20.7,160.0,80.0,6.0,26.43,8.33,21.7,836.0,280.9,5.624,3.26,104.5,70.22,132.3,81.31,702.1],  // MS derived
 ["160 x 80 x 8.0",26.4,160.0,80.0,8.0,33.64,5.0,15.0,1001.0,334.9,5.455,3.155,125.2,83.74,162.9,99.97,882.3],  // MS derived
 ["160 x 80 x 10.0",31.8,160.0,80.0,10.0,40.57,3.0,11.0,1146.0,379.8,5.316,3.06,143.3,94.95,190.9,116.7,1031.0],  // MS derived
 ["180 x 80 x 3.0",11.8,180.0,80.0,3.0,15.01,21.7,55.0,620.8,176.8,6.432,3.432,68.98,44.21,85.82,48.88,444.9],  // MS derived; check availability
 ["180 x 80 x 4.0",15.5,180.0,80.0,4.0,19.75,15.0,40.0,802.1,226.7,6.373,3.388,89.12,56.67,111.8,63.47,578.2],  // MS derived; check availability
 ["180 x 80 x 5.0",19.1,180.0,80.0,5.0,24.36,11.0,31.0,971.0,272.3,6.314,3.343,107.9,68.07,136.5,77.24,704.1],  // MS derived
 ["180 x 80 x 6.0",22.6,180.0,80.0,6.0,28.83,8.33,25.0,1128.0,313.8,6.254,3.299,125.3,78.45,160.0,90.19,822.6],  // MS derived; check availability
 ["180 x 80 x 8.0",28.9,180.0,80.0,8.0,36.84,5.0,17.5,1362.0,376.6,6.079,3.197,151.3,94.15,198.1,111.5,1036.0],  // MS derived
 ["180 x 80 x 10.0",35.0,180.0,80.0,10.0,44.57,3.0,13.0,1570.0,429.1,5.936,3.103,174.5,107.3,233.5,130.7,1214.0],  // MS derived; check availability
 ["180 x 100 x 4.0",16.8,180.0,100.0,4.0,21.35,20.0,40.0,926.0,373.9,6.586,4.185,102.9,74.78,125.9,84.02,853.8],  // MS derived; check availability
 ["180 x 100 x 5.0",20.7,180.0,100.0,5.0,26.36,15.0,31.0,1124.0,451.8,6.531,4.14,124.9,90.35,154.0,102.6,1045.0],  // MS derived
 ["180 x 100 x 6.0",24.5,180.0,100.0,6.0,31.23,11.7,25.0,1310.0,523.8,6.475,4.095,145.5,104.8,180.8,120.2,1227.0],  // MS derived
 ["180 x 100 x 8.0",31.4,180.0,100.0,8.0,40.04,7.5,17.5,1598.0,637.5,6.318,3.99,177.6,127.5,225.6,149.9,1565.0],  // MS derived
 ["180 x 100 x 10.0",38.1,180.0,100.0,10.0,48.57,5.0,13.0,1859.0,736.4,6.188,3.894,206.6,147.3,267.5,177.2,1859.0],  // MS derived
 ["200 x 100 x 4.0",18.0,200.0,100.0,4.0,22.95,20.0,45.0,1200.0,410.8,7.23,4.231,120.0,82.16,148.0,91.7,985.4],  // MS derived
 ["200 x 100 x 5.0",22.3,200.0,100.0,5.0,28.36,15.0,35.0,1459.0,496.9,7.174,4.186,145.9,99.39,181.4,112.1,1206.0],  // MS derived
 ["200 x 100 x 6.0",26.4,200.0,100.0,6.0,33.63,11.7,28.3,1703.0,576.9,7.116,4.142,170.3,115.4,213.3,131.5,1417.0],  // MS derived
 ["200 x 100 x 8.0",33.9,200.0,100.0,8.0,43.24,7.5,20.0,2091.0,705.4,6.954,4.039,209.1,141.1,267.3,164.7,1811.0],  // MS derived
 ["200 x 100 x 10.0",41.3,200.0,100.0,10.0,52.57,5.0,15.0,2444.0,817.7,6.819,3.944,244.4,163.5,318.1,195.2,2154.0],  // MS derived
 ["200 x 120 x 4.0",19.3,200.0,120.0,4.0,24.55,25.0,45.0,1353.0,617.7,7.425,5.016,135.3,102.9,163.7,115.4,1345.0],  // MS derived; check availability
 ["200 x 120 x 5.0",23.8,200.0,120.0,5.0,30.36,19.0,35.0,1649.0,750.1,7.371,4.971,164.9,125.0,200.9,141.4,1652.0],  // MS derived
 ["200 x 120 x 6.0",28.3,200.0,120.0,6.0,36.03,15.0,28.3,1929.0,874.3,7.317,4.926,192.9,145.7,236.5,166.3,1947.0],  // MS derived
 ["200 x 120 x 8.0",36.5,200.0,120.0,8.0,46.44,10.0,20.0,2386.0,1079.0,7.168,4.82,238.6,179.8,298.0,209.5,2507.0],  // MS derived
 ["200 x 120 x 10.0",44.4,200.0,120.0,10.0,56.57,7.0,15.0,2806.0,1262.0,7.043,4.724,280.6,210.4,356.1,249.8,3007.0],  // MS derived
 ["200 x 150 x 4.0",21.2,200.0,150.0,4.0,26.95,32.5,45.0,1584.0,1021.0,7.667,6.155,158.4,136.1,187.2,154.1,1942.0],  // MS derived; check availability
 ["200 x 150 x 5.0",26.2,200.0,150.0,5.0,33.36,25.0,35.0,1935.0,1245.0,7.616,6.109,193.5,166.0,230.1,189.2,2391.0],  // MS derived; check availability
 ["200 x 150 x 6.0",31.1,200.0,150.0,6.0,39.63,20.0,28.3,2268.0,1457.0,7.565,6.063,226.8,194.3,271.5,223.1,2826.0],  // MS derived
 ["200 x 150 x 8.0",40.2,200.0,150.0,8.0,51.24,13.8,20.0,2829.0,1816.0,7.43,5.952,282.9,242.1,344.1,282.8,3665.0],  // MS derived
 ["200 x 150 x 10.0",49.1,200.0,150.0,10.0,62.57,10.0,15.0,3348.0,2143.0,7.315,5.853,334.8,285.8,413.1,339.2,4428.0],  // MS derived
 ["250 x 150 x 5.0",30.1,250.0,150.0,5.0,38.36,25.0,45.0,3304.0,1508.0,9.281,6.27,264.3,201.1,319.8,225.5,3285.0],  // MS derived
 ["250 x 150 x 6.0",35.8,250.0,150.0,6.0,45.63,20.0,36.7,3886.0,1768.0,9.228,6.225,310.8,235.8,378.0,266.3,3886.0],  // MS derived
 ["250 x 150 x 6.3",37.2,250.0,150.0,6.3,47.45,18.8,34.7,4001.0,1825.0,9.183,6.201,320.1,243.3,390.9,275.7,4078.0],  // MS derived; check availability
 ["250 x 150 x 8.0",46.5,250.0,150.0,8.0,59.24,13.8,26.2,4886.0,2219.0,9.081,6.12,390.9,295.9,482.2,339.6,5050.0],  // MS derived
 ["250 x 150 x 10.0",57.0,250.0,150.0,10.0,72.57,10.0,20.0,5825.0,2634.0,8.959,6.025,466.0,351.2,582.0,409.2,6121.0],  // MS derived
 ["250 x 150 x 12.0",66.0,250.0,150.0,12.0,84.06,7.5,15.8,6458.0,2925.0,8.765,5.899,516.6,390.0,658.0,463.3,7088.0],  // MS derived; check availability
 ["250 x 150 x 12.5",68.3,250.0,150.0,12.5,87.04,7.0,15.0,6633.0,3002.0,8.729,5.873,530.6,400.3,678.3,477.5,7315.0],  // MS derived
 ["300 x 100 x 6.0",35.8,300.0,100.0,6.0,45.63,11.7,45.0,4777.0,842.4,10.23,4.296,318.5,168.5,411.4,187.9,2403.0],  // MS derived
 ["300 x 100 x 8.0",46.5,300.0,100.0,8.0,59.24,7.5,32.5,5978.0,1045.0,10.05,4.199,398.5,209.0,523.5,238.3,3080.0],  // MS derived
 ["300 x 100 x 10.0",57.0,300.0,100.0,10.0,72.57,5.0,25.0,7106.0,1224.0,9.896,4.108,473.7,244.9,630.9,285.2,3681.0],  // MS derived
 ["300 x 100 x 12.0",66.0,300.0,100.0,12.0,84.06,3.33,20.0,7808.0,1343.0,9.638,3.997,520.6,268.6,710.3,320.9,4177.0],  // MS derived
 ["300 x 100 x 12.5",68.3,300.0,100.0,12.5,87.04,3.0,19.0,8010.0,1374.0,9.593,3.973,534.0,274.8,731.8,330.2,4292.0],  // MS derived
 ["300 x 200 x 6.0",45.2,300.0,200.0,6.0,57.63,28.3,45.0,7370.0,3962.0,11.31,8.291,491.3,396.2,587.8,446.1,8115.0],  // MS derived
 ["300 x 200 x 6.3",47.1,300.0,200.0,6.3,60.05,26.7,42.6,7624.0,4104.0,11.27,8.267,508.3,410.4,609.9,463.2,8524.0],  // MS derived; check availability
 ["300 x 200 x 8.0",59.1,300.0,200.0,8.0,75.24,20.0,32.5,9389.0,5042.0,11.17,8.186,626.0,504.2,757.1,574.5,10630.0],  // MS derived
 ["300 x 200 x 10.0",72.7,300.0,200.0,10.0,92.57,15.0,25.0,11310.0,6058.0,11.05,8.09,754.2,605.8,920.9,698.1,12990.0],  // MS derived
 ["300 x 200 x 12.0",84.8,300.0,200.0,12.0,108.1,11.7,20.0,12790.0,6854.0,10.88,7.964,852.5,685.4,1056.0,801.2,15240.0],  // MS derived; check availability
 ["300 x 200 x 12.5",88.0,300.0,200.0,12.5,112.0,11.0,19.0,13180.0,7060.0,10.85,7.938,878.6,706.0,1091.0,827.9,15770.0],  // MS derived
 ["400 x 200 x 6.0",54.4,400.0,200.0,6.0,69.63,28.3,61.7,14790.0,5092.0,14.57,8.551,739.5,509.2,906.0,562.5,12070.0],  // MS derived
 ["400 x 200 x 6.3",57.0,400.0,200.0,6.3,72.65,26.7,58.5,15330.0,5286.0,14.53,8.53,766.5,528.6,941.7,585.2,12670.0],  // MS derived; check availability
 ["400 x 200 x 8.0",71.6,400.0,200.0,8.0,91.24,20.0,45.0,18970.0,6517.0,14.42,8.451,948.7,651.7,1173.0,728.1,15820.0],  // MS derived
 ["400 x 200 x 10.0",88.4,400.0,200.0,10.0,112.6,15.0,35.0,23000.0,7864.0,14.3,8.358,1150.0,786.4,1434.0,888.1,19370.0],  // MS derived
 ["400 x 200 x 12.0",104.0,400.0,200.0,12.0,132.1,11.7,28.3,26250.0,8977.0,14.1,8.245,1312.0,897.7,1656.0,1027.0,22780.0],  // MS derived; check availability
 ["400 x 200 x 12.5",108.0,400.0,200.0,12.5,137.0,11.0,27.0,27100.0,9260.0,14.06,8.22,1355.0,926.0,1714.0,1062.0,23590.0],  // MS derived
 ["450 x 250 x 6.0",63.8,450.0,250.0,6.0,81.63,36.7,70.0,22720.0,9245.0,16.68,10.64,1010.0,739.6,1221.0,817.2,20690.0],  // MS derived; check availability
 ["450 x 250 x 6.3",66.9,450.0,250.0,6.3,85.25,34.7,66.4,23610.0,9615.0,16.64,10.62,1049.0,769.2,1271.0,851.5,21730.0],  // MS derived; check availability
 ["450 x 250 x 8.0",84.2,450.0,250.0,8.0,107.2,26.2,51.2,29340.0,11920.0,16.54,10.54,1304.0,953.3,1588.0,1063.0,27220.0],  // MS derived; check availability
 ["450 x 250 x 10.0",104.0,450.0,250.0,10.0,132.6,20.0,40.0,35740.0,14470.0,16.42,10.45,1588.0,1158.0,1948.0,1302.0,33470.0],  // MS derived; check availability
 ["450 x 250 x 12.0",123.0,450.0,250.0,12.0,156.1,15.8,32.5,41140.0,16660.0,16.24,10.33,1828.0,1333.0,2264.0,1515.0,39590.0],  // MS derived; check availability
 ["450 x 250 x 12.5",127.0,450.0,250.0,12.5,162.0,15.0,31.0,42540.0,17220.0,16.2,10.31,1890.0,1377.0,2346.0,1569.0,41060.0],  // MS derived; check availability
 ["500 x 300 x 6.0",73.3,500.0,300.0,6.0,93.63,45.0,78.3,33010.0,15150.0,18.78,12.72,1320.0,1010.0,1581.0,1117.0,32420.0],  // MS derived; check availability
 ["500 x 300 x 6.3",76.8,500.0,300.0,6.3,97.85,42.6,74.4,34350.0,15780.0,18.74,12.7,1374.0,1052.0,1647.0,1165.0,34060.0],  // MS derived; check availability
 ["500 x 300 x 8.0",96.7,500.0,300.0,8.0,123.2,32.5,57.5,42810.0,19620.0,18.64,12.62,1712.0,1308.0,2063.0,1458.0,42770.0],  // MS derived; check availability
 ["500 x 300 x 10.0",120.0,500.0,300.0,10.0,152.6,25.0,45.0,52330.0,23930.0,18.52,12.52,2093.0,1596.0,2537.0,1791.0,52740.0],  // MS derived; check availability
 ["500 x 300 x 12.0",141.0,500.0,300.0,12.0,180.1,20.0,36.7,60600.0,27730.0,18.35,12.41,2424.0,1848.0,2962.0,2093.0,62580.0],  // MS derived; check availability
 ["500 x 300 x 12.5",147.0,500.0,300.0,12.5,187.0,19.0,35.0,62730.0,28690.0,18.31,12.38,2509.0,1912.0,3071.0,2169.0,64950.0]  // MS derived; check availability
].map(a=>({key:a[0],mass:a[1],D:a[2],B:a[3],t:a[4],A:a[5],bT:a[6],dt:a[7],
  Ix:a[8],Iy:a[9],rx:a[10],ry:a[11],Zx:a[12],Zy:a[13],Sx:a[14],Sy:a[15],J:a[16]}));

const RHS_CFmap=Object.fromEntries(RHS_CF.map(s=>[s.key,s]));

const TP385_RHS={"50 x 30 x 3.2":[14.2,6.8],"50 x 30 x 4.0":[16.6,7.77],"50 x 30 x 5.0":[19,8.67],"60 x 40 x 3.2":[30.8,11.7],"60 x 40 x 4.0":[36.7,13.7],"60 x 40 x 5.0":[43,15.7],"60 x 40 x 6.3":[49.5,17.6],"80 x 40 x 3.2":[46.2,16.1],"80 x 40 x 4.0":[55.2,18.9],"80 x 40 x 5.0":[65.1,21.9],"80 x 40 x 6.3":[75.6,24.8],"80 x 40 x 8.0":[85.8,27.4],"90 x 50 x 3.2":[80.9,23.6],"90 x 50 x 4.0":[97.5,28],"90 x 50 x 5.0":[116,32.9],"90 x 50 x 6.3":[138,38.1],"90 x 50 x 8.0":[160,43.2],"100 x 50 x 3.2":[93,26.4],"100 x 50 x 4.0":[113,31.4],"100 x 50 x 5.0":[135,36.9],"100 x 50 x 6.3":[160,42.9],"100 x 50 x 8.0":[186,48.9],"100 x 60 x 3.2":[129,32.4],"100 x 60 x 4.0":[156,38.7],"100 x 60 x 5.0":[188,45.9],"100 x 60 x 6.3":[224,53.8],"100 x 60 x 8.0":[265,62.2],"120 x 60 x 4.0":[201,47.1],"120 x 60 x 5.0":[242,56],"120 x 60 x 6.3":[290,65.9],"120 x 60 x 8.0":[344,76.6],"120 x 80 x 4.0":[330,65],"120 x 80 x 5.0":[401,77.9],"120 x 80 x 6.3":[487,92.9],"120 x 80 x 8.0":[587,110],"120 x 80 x 10.0":[688,126],"150 x 100 x 5.0":[807,127],"150 x 100 x 6.3":[986,153],"150 x 100 x 8.0":[1200,183],"150 x 100 x 10.0":[1430,214],"150 x 100 x 12.5":[1680,246],"160 x 80 x 5.0":[600,106],"160 x 80 x 6.3":[730,127],"160 x 80 x 8.0":[880,151],"160 x 80 x 10.0":[1040,175],"160 x 80 x 12.5":[1200,198],"200 x 100 x 5.0":[1200,172],"200 x 100 x 6.3":[1470,208],"200 x 100 x 8.0":[1800,251],"200 x 100 x 10.0":[2160,295],"200 x 100 x 12.5":[2540,341],"200 x 120 x 6.3":[2030,255],"200 x 120 x 8.0":[2490,310],"200 x 120 x 10.0":[3000,367],"200 x 150 x 8.0":[3640,398],"200 x 150 x 10.0":[4410,475],"250 x 150 x 6.3":[4050,413],"250 x 150 x 8.0":[5020,506],"250 x 150 x 10.0":[6090,605],"250 x 150 x 12.5":[7330,717],"250 x 150 x 16.0":[8870,849],"300 x 100 x 8.0":[3070,387],"300 x 100 x 10.0":[3680,458],"300 x 200 x 6.3":[8480,681],"300 x 200 x 8.0":[10600,840],"300 x 200 x 10.0":[12900,1020],"300 x 200 x 12.5":[15700,1220],"300 x 200 x 16.0":[19300,1470],"400 x 200 x 8.0":[15700,1130],"400 x 200 x 10.0":[19300,1380],"400 x 200 x 12.5":[23400,1660],"400 x 200 x 16.0":[28900,2010],"450 x 250 x 8.0":[27100,1630],"450 x 250 x 10.0":[33300,1990],"450 x 250 x 12.5":[40700,2410],"450 x 250 x 16.0":[50500,2950],"500 x 300 x 10.0":[52400,2700],"500 x 300 x 16.0":[80300,4040]};


/* ==== beam-v03 module 5 ==== */
/* UB section data and P385 torsion properties. */
const UB=[
 // key, mass, D, B, tw, tf, r, d, bT, dt, Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, u, x, J, A, Iw(dm6)
 // Rows tagged "// MS": reference software 2025 UK library import (docs/SECTION_LIBRARY.md). Export values to 4 s.f.;
 // Zx = Ix/(D/2), Zy = Iy/(B/2), rx, ry, bT = (B - tw - 2r)/(2 tf) (EC3 flat outstand, as this file) and dt = d/tw derived.
 ["1016 x 305 x 584",584.0,1056.0,314.0,36.0,64.0,30.0,868.1,1.7,24.1,1246000.0,33400.0,40.9,6.7,23600.0,2130.0,28000.0,3480.0,0.869,18.0,7150.0,744.0,81.2],
 ["1016 x 305 x 494",494.0,1036.0,309.0,31.0,54.0,30.0,868.1,2.02,28.0,1028000.0,26800.0,40.4,6.53,19800.0,1740.0,23400.0,2820.0,0.867,20.9,4400.0,629.0,64.0],
 ["1016 x 305 x 487",486.6,1036.1,308.5,30.0,54.1,30.0,867.9,2.02,28.9,1022000.0,26940.0,40.6,6.59,19730.0,1747.0,23200.0,2800.0,0.871,21.0,4282.0,619.9,63.82],  // MS
 ["1016 x 305 x 438",438.0,1026.0,305.0,26.9,49.0,30.0,868.1,2.23,32.3,910000.0,23400.0,40.4,6.49,17700.0,1540.0,20800.0,2470.0,0.869,23.1,3190.0,556.0,56.0],
 ["1016 x 305 x 437",436.9,1025.9,305.4,26.9,49.0,30.0,867.9,2.23,32.3,910100.0,23660.0,40.4,6.52,17740.0,1549.0,20760.0,2468.0,0.872,23.0,3162.0,556.6,55.5],  // MS
 ["1016 x 305 x 415",415.0,1020.0,304.0,26.0,46.0,30.0,868.1,2.37,33.4,853000.0,21700.0,40.2,6.41,16700.0,1430.0,19600.0,2300.0,0.868,24.3,2700.0,529.0,51.1],
 ["1016 x 305 x 393",392.7,1015.9,303.0,24.4,43.9,30.0,868.1,2.49,35.6,808000.0,20500.0,40.2,6.4,15900.0,1350.0,18500.0,2170.0,0.867,25.5,2330.0,500.0,48.4],
 ["1016 x 305 x 350",350.0,1008.0,302.0,21.1,40.0,30.0,868.1,2.76,41.1,723000.0,18500.0,40.3,6.44,14300.0,1220.0,16600.0,1940.0,0.872,27.8,1720.0,445.0,43.3],
 ["1016 x 305 x 349",349.4,1008.1,302.0,21.1,40.0,30.0,868.1,2.76,41.1,723300.0,18680.0,40.3,6.48,14350.0,1237.0,16590.0,1941.0,0.875,27.9,1690.0,445.1,43.02],  // MS
 ["1016 x 305 x 314",314.3,999.9,300.0,19.1,35.9,30.0,868.1,3.08,45.5,644000.0,16200.0,40.1,6.37,12900.0,1080.0,14800.0,1710.0,0.87,30.7,1260.0,400.0,37.7],
 ["1016 x 305 x 272",272.3,990.1,300.0,16.5,31.0,30.0,868.1,3.6,52.6,554000.0,14000.0,40.0,6.35,11200.0,934.0,12800.0,1470.0,0.871,35.0,835.0,347.0,32.2],
 ["1016 x 305 x 249",248.7,980.1,300.0,16.5,26.0,30.0,868.1,4.3,52.6,481000.0,11800.0,39.0,6.09,9820.0,784.0,11300.0,1240.0,0.86,39.8,582.0,317.0,26.8],
 ["1016 x 305 x 222",222.0,970.3,300.0,16.0,21.1,30.0,868.1,5.31,54.3,408000.0,9550.0,38.0,5.81,8410.0,636.0,9810.0,1020.0,0.85,45.7,390.0,283.0,21.5],
 ["914 x 419 x 388",388.0,921.0,420.5,21.4,36.6,24.1,799.6,4.79,37.4,720000.0,45400.0,38.2,9.59,15600.0,2160.0,17700.0,3340.0,0.885,26.7,1730.0,494.0,88.9],
 ["914 x 419 x 343",343.3,911.8,418.5,19.4,32.0,24.1,799.6,5.48,41.2,626000.0,39200.0,37.8,9.46,13700.0,1870.0,15500.0,2890.0,0.884,30.1,1190.0,437.0,75.8],
 ["914 x 305 x 576",576.0,993.0,322.0,36.1,65.0,19.0,824.8,1.91,22.8,1102000.0,36500.0,38.8,7.06,22200.0,2270.0,26300.0,3660.0,0.874,16.8,7130.0,733.0,77.9],
 ["914 x 305 x 521",521.0,981.0,319.0,33.0,58.9,19.0,824.8,2.11,25.0,982000.0,32100.0,38.5,6.96,20000.0,2020.0,23600.0,3240.0,0.873,18.3,5340.0,664.0,67.7],
 ["914 x 305 x 474",474.0,971.0,316.0,30.0,54.1,19.0,824.8,2.29,27.5,886000.0,28700.0,38.3,6.89,18200.0,1810.0,21400.0,2900.0,0.873,19.8,4100.0,604.0,59.8],
 ["914 x 305 x 425",425.0,961.0,313.0,26.9,49.0,19.0,824.8,2.53,30.7,788000.0,25200.0,38.1,6.82,16400.0,1610.0,19100.0,2560.0,0.874,21.8,3030.0,542.0,52.1],
 ["914 x 305 x 381",381.0,951.0,310.0,24.4,43.9,19.0,824.8,2.82,33.8,697000.0,21900.0,37.9,6.72,14600.0,1410.0,17000.0,2240.0,0.873,24.1,2200.0,486.0,44.8],
 ["914 x 305 x 345",345.0,943.0,308.0,22.1,39.9,19.0,824.8,3.11,37.3,626000.0,19500.0,37.7,6.66,13300.0,1270.0,15400.0,2000.0,0.873,26.3,1650.0,440.0,39.6],
 ["914 x 305 x 313",313.0,932.0,309.0,21.1,34.5,19.0,824.8,3.62,39.1,548000.0,17000.0,37.1,6.54,11800.0,1100.0,13600.0,1750.0,0.867,29.7,1160.0,398.0,34.2],
 ["914 x 305 x 289",289.1,926.6,307.7,19.5,32.0,19.1,824.4,3.91,42.3,504000.0,15600.0,37.0,6.51,10900.0,1010.0,12600.0,1600.0,0.868,31.9,926.0,368.0,31.2],
 ["914 x 305 x 271",271.0,923.0,307.0,18.4,30.0,19.0,824.8,4.18,44.8,472000.0,14500.0,36.9,6.48,10200.0,946.0,11800.0,1490.0,0.867,33.8,769.0,346.0,28.8],
 ["914 x 305 x 253",253.4,918.4,305.5,17.3,27.9,19.1,824.4,4.48,47.7,436000.0,13300.0,36.8,6.42,9500.0,871.0,10900.0,1370.0,0.864,36.2,626.0,323.0,26.4],
 ["914 x 305 x 238",238.0,915.0,305.0,16.5,25.9,19.0,824.8,4.84,50.0,406000.0,12300.0,36.6,6.36,8880.0,806.0,10200.0,1270.0,0.865,38.6,514.0,304.0,24.2],
 ["914 x 305 x 224",224.2,910.4,304.1,15.9,23.9,19.1,824.4,5.23,51.8,376000.0,11200.0,36.3,6.27,8270.0,739.0,9530.0,1160.0,0.86,41.4,422.0,286.0,22.1],
 ["914 x 305 x 201",200.9,903.0,303.3,15.1,20.2,19.1,824.4,6.19,54.6,325000.0,9420.0,35.7,6.07,7200.0,621.0,8350.0,982.0,0.853,46.9,291.0,256.0,18.4],
 ["838 x 292 x 226",226.5,850.9,293.8,16.1,26.8,17.8,761.7,4.52,47.3,340000.0,11400.0,34.3,6.27,7980.0,773.0,9160.0,1210.0,0.87,34.9,514.0,289.0,19.3],
 ["838 x 292 x 194",193.8,840.7,292.4,14.7,21.7,17.8,761.7,5.58,51.8,279000.0,9070.0,33.6,6.06,6640.0,620.0,7640.0,974.0,0.862,41.6,306.0,247.0,15.2],
 ["838 x 292 x 176",175.9,834.9,291.7,14.0,18.8,17.8,761.7,6.44,54.4,246000.0,7800.0,33.1,5.9,5890.0,535.0,6810.0,842.0,0.856,46.5,221.0,224.0,13.0],
 ["762 x 267 x 197",196.8,769.8,268.0,15.6,25.4,16.5,686.0,4.32,44.0,240000.0,8170.0,30.9,5.71,6230.0,610.0,7170.0,958.0,0.869,33.2,404.0,251.0,11.3],
 ["762 x 267 x 173",173.0,762.2,266.7,14.3,21.6,16.5,686.0,5.08,48.0,205000.0,6850.0,30.5,5.58,5390.0,514.0,6200.0,807.0,0.865,38.0,267.0,220.0,9.39],
 ["762 x 267 x 147",146.9,754.0,265.2,12.8,17.5,16.5,686.0,6.27,53.6,169000.0,5460.0,30.0,5.4,4470.0,411.0,5160.0,647.0,0.859,45.2,159.0,187.0,7.4],
 ["762 x 267 x 134",133.9,750.0,264.4,12.0,15.5,16.5,686.0,7.08,57.2,151000.0,4790.0,29.7,5.3,4020.0,362.0,4640.0,570.0,0.853,49.8,119.0,171.0,6.46],
 ["686 x 254 x 170",170.2,692.9,255.8,14.5,23.7,15.2,615.1,4.45,42.4,170000.0,6630.0,28.0,5.53,4920.0,518.0,5630.0,811.0,0.872,31.8,308.0,217.0,7.42],
 ["686 x 254 x 152",152.4,687.5,254.5,13.2,21.0,15.2,615.1,5.02,46.6,150000.0,5780.0,27.8,5.46,4370.0,455.0,5000.0,710.0,0.871,35.4,220.0,194.0,6.42],
 ["686 x 254 x 140",140.1,683.5,253.7,12.4,19.0,15.2,615.1,5.55,49.6,136000.0,5180.0,27.6,5.39,3990.0,409.0,4560.0,638.0,0.87,38.6,169.0,178.0,5.72],
 ["686 x 254 x 125",125.2,677.9,253.0,11.7,16.2,15.2,615.1,6.51,52.6,118000.0,4380.0,27.2,5.24,3480.0,346.0,3990.0,542.0,0.862,43.8,116.0,159.0,4.8],
 ["610 x 305 x 238",238.1,635.8,311.4,18.4,31.4,16.5,540.0,4.14,29.3,209000.0,15800.0,26.3,7.23,6590.0,1020.0,7490.0,1570.0,0.886,21.3,785.0,303.0,14.5],
 ["610 x 305 x 179",179.0,620.2,307.1,14.1,23.6,16.5,540.0,5.51,38.3,153000.0,11400.0,25.9,7.07,4930.0,743.0,5550.0,1140.0,0.885,27.7,340.0,228.0,10.2],
 ["610 x 305 x 149",149.2,612.4,304.8,11.8,19.7,16.5,540.0,6.6,45.8,126000.0,9310.0,25.7,7.0,4110.0,611.0,4590.0,937.0,0.886,32.7,200.0,190.0,8.17],
 ["610 x 229 x 140",139.9,617.2,230.2,13.1,22.1,12.7,547.6,4.34,41.8,112000.0,4510.0,25.0,5.03,3620.0,391.0,4140.0,611.0,0.875,30.5,216.0,178.0,3.99],
 ["610 x 229 x 125",125.1,612.2,229.0,11.9,19.6,12.7,547.6,4.89,46.0,98600.0,3930.0,24.9,4.97,3220.0,343.0,3680.0,535.0,0.875,34.1,154.0,159.0,3.45],
 ["610 x 229 x 113",113.0,607.6,228.2,11.1,17.3,12.7,547.6,5.54,49.3,87300.0,3430.0,24.6,4.88,2870.0,301.0,3280.0,469.0,0.87,38.0,111.0,144.0,2.99],
 ["610 x 229 x 101",101.2,602.6,227.6,10.5,14.8,12.7,547.6,6.48,52.2,75800.0,2910.0,24.2,4.75,2520.0,256.0,2880.0,400.0,0.863,43.1,77.0,129.0,2.52],
 ["610 x 178 x 100",100.3,607.4,179.2,11.3,17.2,12.7,547.6,4.14,48.5,72500.0,1660.0,23.8,3.6,2390.0,185.0,2790.0,296.0,0.855,38.7,95.0,128.0,1.44],
 ["610 x 178 x 92",92.2,603.0,178.8,10.9,15.0,12.7,547.6,4.75,50.2,64600.0,1440.0,23.4,3.5,2140.0,161.0,2510.0,258.0,0.85,42.6,71.0,117.0,1.24],
 ["610 x 178 x 82",81.8,598.6,177.9,10.0,12.8,12.7,547.6,5.57,54.8,55900.0,1210.0,23.2,3.4,1870.0,136.0,2190.0,218.0,0.843,48.4,48.8,104.0,1.04],
 ["533 x 312 x 273",273.3,577.1,320.2,21.1,37.6,12.7,476.5,3.64,22.6,199000.0,20600.0,23.9,7.69,6890.0,1290.0,7870.0,1990.0,0.891,15.9,1290.0,348.0,15.0],
 ["533 x 312 x 272",273.2,577.1,320.2,21.1,37.6,12.7,476.5,3.64,22.6,198600.0,20620.0,23.9,7.7,6883.0,1288.0,7858.0,1985.0,0.891,15.9,1288.0,348.1,14.97],  // MS; EN 10365 designation, identical geometry to "533 x 312 x 273"
 ["533 x 312 x 219",218.8,560.3,317.4,18.3,29.2,12.7,476.5,4.69,26.0,151000.0,15600.0,23.3,7.48,5400.0,982.0,6120.0,1510.0,0.884,19.8,642.0,279.0,11.0],
 ["533 x 312 x 182",181.5,550.7,314.5,15.2,24.4,12.7,476.5,5.61,31.3,123000.0,12700.0,23.1,7.4,4480.0,806.0,5040.0,1240.0,0.887,23.4,373.0,231.0,8.77],
 ["533 x 312 x 151",150.6,542.5,312.0,12.7,20.3,12.7,476.5,6.75,37.5,101000.0,10300.0,22.9,7.32,3710.0,659.0,4150.0,1010.0,0.886,27.8,216.0,192.0,7.01],
 ["533 x 312 x 150",150.6,542.5,312.0,12.7,20.3,12.7,476.5,6.75,37.5,100600.0,10290.0,22.9,7.32,3709.0,659.6,4142.0,1010.0,0.885,27.8,216.2,191.8,7.005],  // MS; EN 10365 designation, identical geometry to "533 x 312 x 151"
 ["533 x 210 x 138",138.3,549.1,213.9,14.7,23.6,12.7,476.5,3.68,32.4,86100.0,3860.0,22.1,4.68,3140.0,361.0,3610.0,568.0,0.873,25.0,250.0,176.0,2.67],
 ["533 x 210 x 122",122.0,544.5,211.9,12.7,21.3,12.7,476.5,4.08,37.5,76000.0,3390.0,22.1,4.67,2790.0,320.0,3200.0,500.0,0.878,27.6,178.0,155.0,2.32],
 ["533 x 210 x 109",109.0,539.5,210.8,11.6,18.8,12.7,476.5,4.62,41.1,66800.0,2940.0,21.9,4.6,2480.0,279.0,2830.0,436.0,0.875,30.9,126.0,139.0,1.99],
 ["533 x 210 x 101",101.0,536.7,210.0,10.8,17.4,12.7,476.5,4.99,44.1,61500.0,2690.0,21.9,4.57,2290.0,256.0,2610.0,399.0,0.873,33.2,101.0,129.0,1.81],
 ["533 x 210 x 92",92.1,533.1,209.3,10.1,15.6,12.7,476.5,5.57,47.2,55200.0,2390.0,21.7,4.51,2070.0,228.0,2360.0,355.0,0.873,36.4,75.7,117.0,1.6],
 ["533 x 210 x 82",82.2,528.3,208.8,9.6,13.2,12.7,476.5,6.58,49.6,47500.0,2010.0,21.3,4.38,1800.0,192.0,2060.0,300.0,0.864,41.5,51.5,105.0,1.33],
 ["533 x 165 x 85",84.8,534.9,166.5,10.3,16.5,12.7,476.5,3.96,46.3,48500.0,1270.0,21.2,3.44,1820.0,153.0,2100.0,243.0,0.859,35.5,73.8,108.0,0.857],
 ["533 x 165 x 75",74.7,529.1,165.9,9.7,13.6,12.7,476.5,4.81,49.1,41100.0,1040.0,20.8,3.3,1550.0,125.0,1810.0,200.0,0.853,41.1,47.9,95.2,0.691],
 ["533 x 165 x 74",74.7,529.1,165.9,9.7,13.6,12.7,476.5,4.81,49.1,41060.0,1040.0,20.8,3.31,1552.0,125.4,1808.0,200.0,0.854,41.0,47.94,95.19,0.6876],  // MS; EN 10365 designation, identical geometry to "533 x 165 x 75"
 ["533 x 165 x 66",65.7,524.7,165.1,8.9,11.4,12.7,476.5,5.74,53.5,35000.0,859.0,20.5,3.2,1340.0,104.0,1560.0,166.0,0.847,47.0,32.0,83.7,0.566],
 ["457 x 191 x 161",161.4,492.0,199.4,18.0,32.0,10.2,407.6,2.52,22.6,79800.0,4250.0,19.7,4.55,3240.0,426.0,3780.0,672.0,0.881,16.5,515.0,206.0,2.25],
 ["457 x 191 x 133",133.3,480.6,196.7,15.3,26.3,10.2,407.6,3.06,26.6,63800.0,3350.0,19.4,4.44,2660.0,341.0,3070.0,535.0,0.88,19.6,292.0,170.0,1.73],
 ["457 x 191 x 106",105.8,469.2,194.0,12.6,20.6,10.2,407.6,3.91,32.3,48900.0,2510.0,19.0,4.32,2080.0,259.0,2390.0,405.0,0.876,24.5,146.0,135.0,1.27],
 ["457 x 191 x 98",98.3,467.2,192.8,11.4,19.6,10.2,407.6,4.11,35.8,45700.0,2350.0,19.1,4.33,1960.0,243.0,2230.0,379.0,0.881,25.8,121.0,125.0,1.18],
 ["457 x 191 x 89",89.3,463.4,191.9,10.5,17.7,10.2,407.6,4.55,38.8,41000.0,2090.0,19.0,4.29,1770.0,218.0,2010.0,338.0,0.877,28.3,90.7,114.0,1.04],
 ["457 x 191 x 82",82.0,460.0,191.3,9.9,16.0,10.2,407.6,5.03,41.2,37100.0,1870.0,18.8,4.23,1610.0,196.0,1830.0,304.0,0.879,30.8,69.2,104.0,0.922],
 ["457 x 191 x 74",74.3,457.0,190.4,9.0,14.5,10.2,407.6,5.55,45.3,33300.0,1670.0,18.8,4.2,1460.0,176.0,1650.0,272.0,0.876,33.8,51.8,94.6,0.818],
 ["457 x 191 x 67",67.1,453.4,189.9,8.5,12.7,10.2,407.6,6.34,48.0,29400.0,1450.0,18.5,4.12,1300.0,153.0,1470.0,237.0,0.872,37.9,37.1,85.5,0.705],
 ["457 x 152 x 82",82.1,465.8,155.3,10.5,18.9,10.2,407.6,3.29,38.8,36600.0,1180.0,18.7,3.37,1570.0,153.0,1810.0,240.0,0.87,27.5,89.2,105.0,0.591],
 ["457 x 152 x 74",74.2,462.0,154.4,9.6,17.0,10.2,407.6,3.66,42.5,32700.0,1050.0,18.6,3.33,1410.0,136.0,1630.0,213.0,0.874,30.1,65.9,94.5,0.518],
 ["457 x 152 x 67",67.2,458.0,153.8,9.0,15.0,10.2,407.6,4.15,45.3,28900.0,913.0,18.4,3.27,1260.0,119.0,1450.0,187.0,0.867,33.6,47.7,85.6,0.448],
 ["457 x 152 x 60",59.8,454.6,152.9,8.1,13.3,10.2,407.6,4.68,50.3,25500.0,795.0,18.3,3.23,1120.0,104.0,1290.0,163.0,0.869,37.5,33.8,76.2,0.387],
 ["457 x 152 x 52",52.3,449.8,152.4,7.6,10.9,10.2,407.6,5.71,53.6,21400.0,645.0,17.9,3.11,950.0,84.6,1100.0,133.0,0.861,43.8,21.4,66.6,0.311],
 ["406 x 178 x 85",85.3,417.2,181.9,10.9,18.2,10.2,360.4,4.14,33.1,31700.0,1830.0,17.1,4.11,1520.0,201.0,1730.0,313.0,0.879,24.4,93.0,109.0,0.728],
 ["406 x 178 x 74",74.2,412.8,179.5,9.5,16.0,10.2,360.4,4.68,37.9,27300.0,1550.0,17.0,4.04,1320.0,172.0,1500.0,267.0,0.882,27.5,62.8,94.5,0.608],
 ["406 x 178 x 67",67.1,409.4,178.8,8.8,14.3,10.2,360.4,5.23,41.0,24300.0,1360.0,16.9,3.99,1190.0,153.0,1350.0,237.0,0.88,30.5,46.1,85.5,0.533],
 ["406 x 178 x 60",60.1,406.4,177.9,7.9,12.8,10.2,360.4,5.84,45.6,21600.0,1200.0,16.8,3.97,1060.0,135.0,1200.0,209.0,0.88,33.8,33.3,76.5,0.466],
 ["406 x 178 x 54",54.1,402.6,177.7,7.7,10.9,10.2,360.4,6.86,46.8,18700.0,1020.0,16.5,3.85,930.0,115.0,1050.0,178.0,0.869,38.3,23.1,69.0,0.392],
 ["406 x 140 x 53",53.3,406.6,143.3,7.9,12.9,10.2,360.4,4.46,45.6,18300.0,635.0,16.4,3.06,899.0,88.6,1030.0,139.0,0.87,34.1,29.0,67.9,0.246],
 ["406 x 140 x 46",46.0,403.2,142.2,6.8,11.2,10.2,360.4,5.13,53.0,15700.0,538.0,16.4,3.03,778.0,75.7,888.0,118.0,0.871,39.0,19.0,58.6,0.207],
 ["406 x 140 x 39",39.0,398.0,141.8,6.4,8.6,10.2,360.4,6.69,56.3,12500.0,410.0,15.9,2.87,629.0,57.8,724.0,90.8,0.858,47.4,10.7,49.7,0.155],
 ["356 x 171 x 67",67.1,363.4,173.2,9.1,15.7,10.2,311.6,4.58,34.2,19500.0,1360.0,15.1,3.99,1070.0,157.0,1210.0,243.0,0.886,24.4,55.7,85.5,0.412],
 ["356 x 171 x 57",57.0,358.0,172.2,8.1,13.0,10.2,311.6,5.53,38.5,16000.0,1110.0,14.9,3.91,896.0,129.0,1010.0,199.0,0.882,28.8,33.4,72.6,0.33],
 ["356 x 171 x 51",51.0,355.0,171.5,7.4,11.5,10.2,311.6,6.25,42.1,14100.0,968.0,14.8,3.86,796.0,113.0,896.0,174.0,0.88,32.1,23.8,64.9,0.286],
 ["356 x 171 x 45",45.0,351.4,171.1,7.0,9.7,10.2,311.6,7.41,44.5,12100.0,811.0,14.5,3.76,687.0,94.8,775.0,147.0,0.874,36.8,15.8,57.3,0.237],
 ["356 x 127 x 39",39.1,353.4,126.0,6.6,10.7,10.2,311.6,4.63,47.2,10200.0,358.0,14.3,2.68,576.0,56.8,659.0,89.0,0.871,35.2,15.1,49.8,0.105],
 ["356 x 127 x 33",33.1,349.0,125.4,6.0,8.5,10.2,311.6,5.82,51.9,8250.0,280.0,14.0,2.58,473.0,44.7,543.0,70.2,0.863,42.1,8.79,42.1,0.081],
 ["305 x 165 x 54",54.0,310.4,166.9,7.9,13.7,8.9,265.2,5.15,33.6,11700.0,1060.0,13.0,3.93,754.0,127.0,846.0,196.0,0.888,23.6,34.8,68.8,0.234],
 ["305 x 165 x 46",46.1,306.6,165.7,6.7,11.8,8.9,265.2,5.98,39.6,9900.0,896.0,13.0,3.9,646.0,108.0,720.0,166.0,0.89,27.1,22.2,58.7,0.195],
 ["305 x 165 x 40",40.3,303.4,165.0,6.0,10.2,8.9,265.2,6.92,44.2,8500.0,764.0,12.9,3.86,560.0,92.6,623.0,142.0,0.889,31.0,14.7,51.3,0.164],
 ["305 x 127 x 48",48.1,311.0,125.3,9.0,14.0,8.9,265.2,3.52,29.5,9570.0,461.0,12.5,2.74,616.0,73.6,711.0,116.0,0.873,23.3,31.8,61.2,0.102],
 ["305 x 127 x 42",41.9,307.2,124.3,8.0,12.1,8.9,265.2,4.07,33.2,8200.0,389.0,12.4,2.7,534.0,62.6,614.0,98.4,0.872,26.5,21.1,53.4,0.0846],
 ["305 x 127 x 37",37.0,304.4,123.4,7.1,10.7,8.9,265.2,4.6,37.4,7170.0,336.0,12.3,2.67,471.0,54.5,539.0,85.4,0.871,29.7,14.8,47.2,0.0725],
 ["305 x 102 x 33",32.8,312.7,102.4,6.6,10.8,7.6,275.9,3.73,41.8,6500.0,194.0,12.5,2.15,416.0,37.9,481.0,60.0,0.867,31.6,12.2,41.8,0.0442],
 ["305 x 102 x 28",28.2,308.7,101.8,6.0,8.8,7.6,275.9,4.58,46.0,5370.0,155.0,12.2,2.08,348.0,30.5,403.0,48.4,0.859,37.4,7.4,35.9,0.0349],
 ["305 x 102 x 25",24.8,305.1,101.6,5.8,7.0,7.6,275.9,5.76,47.6,4460.0,123.0,11.9,1.97,292.0,24.2,342.0,38.8,0.849,43.1,4.77,31.6,0.027],
 ["254 x 146 x 43",43.0,259.6,147.3,7.2,12.7,7.6,219.0,4.92,30.4,6540.0,677.0,10.9,3.52,504.0,92.0,566.0,141.0,0.89,21.1,23.9,54.8,0.103],
 ["254 x 146 x 37",37.0,256.0,146.4,6.3,10.9,7.6,219.0,5.73,34.8,5540.0,571.0,10.8,3.48,433.0,78.0,483.0,119.0,0.889,24.3,15.3,47.2,0.0857],
 ["254 x 146 x 31",31.1,251.4,146.1,6.0,8.6,7.6,219.0,7.26,36.5,4410.0,448.0,10.5,3.36,351.0,61.3,393.0,94.1,0.879,29.6,8.55,39.7,0.066],
 ["254 x 102 x 28",28.3,260.4,102.2,6.3,10.0,7.6,225.2,4.04,35.7,4000.0,179.0,10.5,2.22,308.0,34.9,353.0,54.8,0.874,27.5,9.57,36.1,0.028],
 ["254 x 102 x 25",25.2,257.2,101.9,6.0,8.4,7.6,225.2,4.8,37.5,3410.0,149.0,10.3,2.15,266.0,29.2,306.0,46.0,0.868,31.4,6.42,32.0,0.023],
 ["254 x 102 x 22",22.0,254.0,101.6,5.7,6.8,7.6,225.2,5.93,39.5,2840.0,119.0,10.1,2.06,224.0,23.5,259.0,37.3,0.856,36.3,4.15,28.0,0.0182],
 ["203 x 133 x 30",30.0,206.8,133.9,6.4,9.6,7.6,172.4,5.85,26.9,2900.0,385.0,8.71,3.17,280.0,57.5,314.0,88.2,0.881,21.5,10.3,38.2,0.0374],
 ["203 x 133 x 25",25.1,203.2,133.2,5.7,7.8,7.6,172.4,7.2,30.2,2340.0,308.0,8.56,3.1,230.0,46.2,258.0,70.9,0.877,25.6,5.96,32.0,0.0294],
 ["203 x 102 x 23",23.1,203.2,101.8,5.4,9.3,7.6,169.4,4.37,31.4,2100.0,164.0,8.46,2.36,207.0,32.2,234.0,49.7,0.888,22.4,7.02,29.4,0.0154],
 ["178 x 102 x 19",19.0,177.8,101.2,4.8,7.9,7.6,146.8,5.14,30.6,1360.0,137.0,7.48,2.37,153.0,27.0,171.0,41.6,0.886,22.6,4.41,24.3,0.0099],
 ["152 x 89 x 16",16.0,152.4,88.7,4.5,7.7,7.6,121.8,4.48,27.1,834.0,89.8,6.41,2.1,109.0,20.2,123.0,31.2,0.889,19.5,3.56,20.3,0.0047],
 ["127 x 76 x 13",13.0,127.0,76.0,4.0,7.6,7.6,96.6,3.74,24.2,473.0,55.7,5.35,1.84,74.6,14.7,84.2,22.6,0.894,16.3,2.85,16.5,0.002]
].map(a=>({key:a[0],mass:a[1],D:a[2],B:a[3],tw:a[4],tf:a[5],r:a[6],d:a[7],bT:a[8],dt:a[9],
  Ix:a[10],Iy:a[11],rx:a[12],ry:a[13],Zx:a[14],Zy:a[15],Sx:a[16],Sy:a[17],u:a[18],x:a[19],J:a[20],A:a[21],Iw:a[22]}));

const UBmap=Object.fromEntries(UB.map(s=>[s.key,s]));

const TP385_UB={"1016 x 305 x 438":[3180,2.13,55.9,746,27900],"1016 x 305 x 393":[2330,2.32,48.4,736,24500],"1016 x 305 x 350":[1720,2.56,43.3,731,22100],"1016 x 305 x 314":[1260,2.78,37.6,723,19500],"1016 x 305 x 272":[835,3.16,32.2,719,16700],"1016 x 305 x 249":[582,3.46,26.9,716,14000],"1016 x 305 x 222":[390,3.78,21.5,712,11300],"914 x 419 x 388":[1730,3.64,88.8,930,35800],"914 x 419 x 343":[1190,4.06,75.9,920,30800],"914 x 305 x 289":[926,2.96,31.2,688,16900],"914 x 305 x 253":[626,3.3,26.4,680,14500],"914 x 305 x 224":[422,3.68,22,674,12200],"914 x 305 x 201":[291,4.04,18.4,669,10300],"838 x 292 x 226":[514,3.13,19.4,605,11900],"838 x 292 x 194":[306,3.59,15.2,599,9500],"838 x 292 x 176":[221,3.9,13,595,8160],"762 x 267 x 197":[404,2.69,11.3,499,8490],"762 x 267 x 173":[267,3.02,9.39,494,7110],"762 x 267 x 147":[159,3.47,7.4,488,5670],"762 x 267 x 134":[119,3.75,6.46,486,4970],"686 x 254 x 170":[308,2.5,7.42,428,6490],"686 x 254 x 152":[220,2.75,6.42,424,5670],"686 x 254 x 140":[169,2.96,5.72,421,5080],"686 x 254 x 125":[116,3.27,4.79,419,4290],"610 x 305 x 238":[785,2.18,14.4,471,11500],"610 x 305 x 179":[340,2.78,10.1,458,8300],"610 x 305 x 149":[200,3.26,8.18,452,6780],"610 x 229 x 140":[216,2.19,3.99,342,4360],"610 x 229 x 125":[154,2.41,3.45,339,3810],"610 x 229 x 113":[111,2.64,2.99,337,3320],"610 x 229 x 101":[77,2.91,2.51,334,2820],"610 x 178 x 100":[95,1.99,1.45,264,2040],"610 x 178 x 92":[71,2.13,1.24,263,1760],"610 x 178 x 82":[48.8,2.35,1.04,261,1480],"533 x 312 x 273":[1290,1.74,15,432,13000],"533 x 312 x 219":[642,2.11,11,422,9770],"533 x 312 x 182":[373,2.47,8.79,414,7950],"533 x 312 x 151":[216,2.9,7.03,408,6460],"533 x 210 x 138":[250,1.66,2.67,281,3550],"533 x 210 x 122":[178,1.84,2.32,277,3130],"533 x 210 x 109":[126,2.02,1.99,274,2720],"533 x 210 x 101":[101,2.16,1.81,273,2490],"533 x 210 x 92":[76,2.34,1.6,271,2210],"533 x 210 x 82":[52,2.59,1.33,269,1850],"533 x 165 x 85":[73.8,1.73,0.85,215,1480],"533 x 165 x 75":[47.9,1.93,0.69,214,1210],"533 x 165 x 66":[32,2.14,0.57,212,997],"457 x 191 x 161":[515,1.06,2.25,229,3660],"457 x 191 x 133":[292,1.24,1.73,223,2890],"457 x 191 x 106":[146,1.5,1.27,218,2170],"457 x 191 x 98":[121,1.59,1.18,216,2040],"457 x 191 x 89":[90.7,1.72,1.04,214,1820],"457 x 191 x 82":[69.2,1.86,0.922,212,1620],"457 x 191 x 74":[51.8,2.02,0.818,211,1450],"457 x 191 x 67":[37.1,2.22,0.705,209,1260],"457 x 152 x 82":[89.2,1.31,0.592,174,1270],"457 x 152 x 74":[65.9,1.43,0.518,172,1130],"457 x 152 x 67":[47.7,1.56,0.448,170,982],"457 x 152 x 60":[33.8,1.72,0.387,169,858],"457 x 152 x 52":[21.4,1.94,0.311,167,694],"406 x 178 x 85":[93,1.43,0.728,181,1500],"406 x 178 x 74":[62.8,1.58,0.608,178,1280],"406 x 178 x 67":[46.1,1.73,0.533,177,1130],"406 x 178 x 60":[33.3,1.9,0.466,175,997],"406 x 178 x 54":[23.1,2.1,0.392,174,843],"406 x 140 x 53":[29,1.48,0.246,141,652],"406 x 140 x 46":[19,1.68,0.207,139,555],"406 x 140 x 39":[10.7,1.94,0.155,138,421],"356 x 171 x 67":[55.7,1.38,0.412,151,1020],"356 x 171 x 57":[33.4,1.6,0.33,149,831],"356 x 171 x 51":[23.8,1.76,0.286,147,726],"356 x 171 x 45":[15.8,1.97,0.237,146,606],"356 x 127 x 39":[15.1,1.34,0.105,108,364],"356 x 127 x 33":[8.79,1.55,0.081,107,284],"305 x 165 x 54":[34.8,1.32,0.234,124,708],"305 x 165 x 46":[22.2,1.51,0.195,122,597],"305 x 165 x 40":[14.7,1.7,0.164,121,509],"305 x 127 x 48":[31.8,0.91,0.102,93,408],"305 x 127 x 42":[21.1,1.02,0.085,91.7,345],"305 x 127 x 37":[14.8,1.13,0.072,90.6,299],"305 x 102 x 33":[12.2,0.97,0.0442,77.3,214],"305 x 102 x 28":[7.4,1.11,0.0349,76.3,171],"305 x 102 x 25":[4.77,1.22,0.0273,75.7,135],"254 x 146 x 43":[23.9,1.06,0.103,90.9,425],"254 x 146 x 37":[15.3,1.2,0.0858,89.7,358],"254 x 146 x 31":[8.55,1.41,0.066,88.7,279],"254 x 102 x 28":[9.57,0.87,0.0281,64,163],"254 x 102 x 25":[6.42,0.97,0.0231,63.4,136],"254 x 102 x 22":[4.15,1.07,0.0182,62.8,108],"203 x 133 x 30":[10.3,0.97,0.0374,66,212],"203 x 133 x 25":[5.96,1.13,0.0294,65.1,169],"203 x 102 x 23":[7.02,0.75,0.0154,49.3,117],"178 x 102 x 19":[4.41,0.76,0.0099,43,85.9],"152 x 89 x 16":[3.56,0.59,0.0047,32.1,54.8],"127 x 76 x 13":[2.85,0.43,0.002,22.7,32.8],
 // Rows below are not in SCI P385 Appendix A. Derived: IT = J and Iw from the table row above,
 // a = sqrt(E*Iw/(G*IT)) with E = 210000, G = 81000 N/mm2 (m), Wn0 = (D - tf)*B/4 (cm2), Sw1 = (D - tf)*B^2*tf/16 (cm4).
 // Formulas reproduce every tabulated P385 row of this family within 0.5 % (Wn0, Sw1) and 0.9 % (a); see docs/SECTION_LIBRARY.md.
 "1016 x 305 x 584":[7150,1.72,81.2,779,39100],  // derived from table row
 "1016 x 305 x 494":[4400,1.94,64,759,31600],  // derived from table row
 "1016 x 305 x 415":[2700,2.22,51.1,740,25900],  // derived from table row
 "914 x 305 x 576":[7130,1.68,77.9,747,39100],  // derived from table row
 "914 x 305 x 521":[5340,1.81,67.7,735,34500],  // derived from table row
 "914 x 305 x 474":[4100,1.94,59.8,724,31000],  // derived from table row
 "914 x 305 x 425":[3030,2.11,52.1,714,27400],  // derived from table row
 "914 x 305 x 381":[2200,2.3,44.8,703,23900],  // derived from table row
 "914 x 305 x 345":[1650,2.49,39.6,695,21400],  // derived from table row
 "914 x 305 x 313":[1160,2.76,34.2,693,18500],  // derived from table row
 "914 x 305 x 271":[769,3.12,28.8,685,15800],  // derived from table row
 "914 x 305 x 238":[514,3.49,24.2,678,13400],  // derived from table row
 "1016 x 305 x 487":[4282,1.97,63.82,757,31600],  // MS, derived
 "1016 x 305 x 437":[3162,2.13,55.5,746,27900],  // MS, derived
 "1016 x 305 x 349":[1690,2.57,43.02,731,22100],  // MS, derived
 "533 x 312 x 272":[1288,1.74,14.97,432,13000],  // MS, derived
 "533 x 312 x 150":[216.2,2.9,7.005,407,6450],  // MS, derived
 "533 x 165 x 74":[47.94,1.93,0.6876,214,1210]  // MS, derived
};


/* ==== beam-v03 module 6 ==== */
/* UC section data and P385 torsion properties. */
const UC=[
 // key, mass, D, B, tw, tf, r, d, bT, dt, Ix, Iy, rx, ry, Zx, Zy, Sx, Sy, u, x, J, A, Iw(dm6)
 // Rows tagged "// MS": reference software 2025 UK library import (docs/SECTION_LIBRARY.md). Export values to 4 s.f.;
 // Zx = Ix/(D/2), Zy = Iy/(B/2), rx, ry, bT = (B - tw - 2r)/(2 tf) (EC3 flat outstand, as this file) and dt = d/tw derived.
 ["356 x 406 x 1299",1299.0,600.0,476.0,100.0,140.0,15.4,290.0,1.23,2.9,755000.0,254000.0,21.4,12.4,25200.0,10700.0,33200.0,16700.0,0.846,3.42,94500.0,1655.0,133.1],
 ["356 x 406 x 1202",1202.0,580.0,471.0,95.0,130.0,15.4,290.0,1.33,3.05,664000.0,229000.0,20.8,12.2,22900.0,9710.0,30000.0,15200.0,0.842,3.59,76300.0,1531.0,114.6],
 ["356 x 406 x 1086",1086.0,569.0,454.0,78.0,125.0,15.0,290.0,1.38,3.72,596000.0,196000.0,20.7,11.9,20900.0,8640.0,27200.0,13400.0,0.852,3.79,60500.0,1386.0,96.1],
 ["356 x 406 x 990",990.0,550.0,448.0,71.9,115.0,15.0,290.0,1.5,4.03,519000.0,173000.0,20.3,11.7,18900.0,7740.0,24300.0,12000.0,0.851,4.02,46900.0,1262.0,81.5],
 ["356 x 406 x 900",900.0,531.0,442.0,65.9,106.0,15.0,290.0,1.63,4.4,450000.0,153000.0,19.8,11.6,17000.0,6940.0,21600.0,10700.0,0.849,4.26,36400.0,1149.0,68.9],
 ["356 x 406 x 818",818.0,514.0,437.0,60.5,97.0,15.0,290.0,1.79,4.79,392000.0,136000.0,19.4,11.4,15300.0,6200.0,19300.0,9560.0,0.847,4.55,27800.0,1043.0,58.6],
 ["356 x 406 x 744",744.0,498.0,432.0,55.6,88.9,15.0,290.0,1.95,5.22,342000.0,120000.0,19.0,11.3,13700.0,5550.0,17200.0,8550.0,0.845,4.86,21400.0,948.0,50.0],
 ["356 x 406 x 677",677.0,483.0,428.0,51.2,81.5,15.0,290.0,2.13,5.66,300000.0,107000.0,18.6,11.1,12400.0,4990.0,15400.0,7680.0,0.844,5.19,16400.0,863.0,42.9],
 ["356 x 406 x 634",633.9,474.6,424.0,47.6,77.0,15.2,290.2,2.25,6.1,275000.0,98100.0,18.4,11.0,11600.0,4630.0,14200.0,7110.0,0.842,5.46,13700.0,808.0,38.8],
 ["356 x 406 x 592",592.0,465.0,421.0,45.0,72.3,15.0,290.0,2.39,6.44,250000.0,90200.0,18.2,10.9,10800.0,4280.0,13100.0,6570.0,0.843,5.72,11400.0,755.0,34.7],
 ["356 x 406 x 551",551.0,455.6,418.5,42.1,67.5,15.2,290.2,2.56,6.89,227000.0,82700.0,18.0,10.9,9960.0,3950.0,12100.0,6060.0,0.842,6.05,9240.0,702.0,31.1],
 ["356 x 406 x 509",509.0,446.0,416.0,39.1,62.7,15.0,290.0,2.77,7.42,204000.0,75400.0,17.8,10.8,9170.0,3620.0,11000.0,5550.0,0.84,6.42,7390.0,649.0,27.6],
 ["356 x 406 x 477",477.0,427.0,424.4,48.0,53.2,15.2,290.2,3.25,6.05,172500.0,68100.0,16.9,10.6,8080.0,3209.0,9704.0,4981.0,0.816,6.89,5700.0,607.4,23.68],  // MS
 ["356 x 406 x 467",467.0,436.6,412.2,35.8,58.0,15.2,290.2,2.98,8.11,183000.0,67800.0,17.5,10.7,8380.0,3290.0,10000.0,5030.0,0.839,6.85,5810.0,595.0,24.3],
 ["356 x 406 x 393",393.0,419.0,407.0,30.6,49.2,15.2,290.2,3.52,9.48,147000.0,55400.0,17.1,10.5,7000.0,2720.0,8220.0,4150.0,0.837,7.85,3550.0,501.0,18.9],
 ["356 x 406 x 340",339.9,406.4,403.0,26.6,42.9,15.2,290.2,4.03,10.9,123000.0,46900.0,16.8,10.4,6030.0,2330.0,7000.0,3540.0,0.836,8.85,2340.0,433.0,15.5],
 ["356 x 406 x 287",287.1,393.6,399.0,22.6,36.5,15.2,290.2,4.74,12.8,99900.0,38700.0,16.5,10.3,5070.0,1940.0,5810.0,2950.0,0.835,10.2,1440.0,366.0,12.3],
 ["356 x 406 x 235",235.1,381.0,394.8,18.4,30.2,15.2,290.2,5.73,15.8,79100.0,31000.0,16.3,10.2,4150.0,1570.0,4690.0,2380.0,0.835,12.0,812.0,299.0,9.54],
 ["356 x 368 x 202",201.9,374.6,374.7,16.5,27.0,15.2,290.2,6.07,17.6,66300.0,23700.0,16.1,9.6,3540.0,1260.0,3970.0,1920.0,0.844,13.3,558.0,257.0,7.16],
 ["356 x 368 x 177",177.0,368.2,372.6,14.4,23.8,15.2,290.2,6.89,20.2,57100.0,20500.0,15.9,9.54,3100.0,1100.0,3460.0,1670.0,0.843,15.0,381.0,226.0,6.09],
 ["356 x 368 x 153",152.9,362.0,370.5,12.3,20.7,15.2,290.2,7.92,23.6,48600.0,17600.0,15.8,9.49,2680.0,948.0,2960.0,1430.0,0.843,17.0,251.0,195.0,5.11],
 ["356 x 368 x 129",129.0,355.6,368.6,10.4,17.5,15.2,290.2,9.37,27.9,40200.0,14600.0,15.6,9.43,2260.0,793.0,2480.0,1200.0,0.845,19.8,153.0,164.0,4.18],
 ["305 x 305 x 283",282.9,365.3,322.2,26.8,44.1,15.2,246.7,3.0,9.21,78900.0,24600.0,14.8,8.27,4320.0,1530.0,5110.0,2340.0,0.856,7.65,2030.0,360.0,6.35],
 ["305 x 305 x 240",240.0,352.5,318.4,23.0,37.7,15.2,246.7,3.51,10.7,64200.0,20300.0,14.5,8.15,3640.0,1280.0,4250.0,1950.0,0.854,8.74,1270.0,306.0,5.03],
 ["305 x 305 x 198",198.1,339.9,314.5,19.1,31.4,15.2,246.7,4.22,12.9,50900.0,16300.0,14.2,8.04,3000.0,1040.0,3440.0,1580.0,0.854,10.2,734.0,252.0,3.88],
 ["305 x 305 x 158",158.1,327.1,311.2,15.8,25.0,15.2,246.7,5.3,15.6,38700.0,12600.0,13.9,7.9,2370.0,808.0,2680.0,1230.0,0.852,12.4,378.0,201.0,2.87],
 ["305 x 305 x 137",136.9,320.5,309.2,13.8,21.7,15.2,246.7,6.11,17.9,32800.0,10700.0,13.7,7.83,2050.0,692.0,2300.0,1050.0,0.852,14.1,249.0,174.0,2.39],
 ["305 x 305 x 118",117.9,314.5,307.4,12.0,18.7,15.2,246.7,7.09,20.6,27700.0,9060.0,13.6,7.77,1760.0,589.0,1960.0,895.0,0.852,16.1,161.0,150.0,1.98],
 ["305 x 305 x 97",96.9,307.9,305.3,9.9,15.4,15.2,246.7,8.6,24.9,22200.0,7310.0,13.4,7.69,1450.0,479.0,1590.0,726.0,0.851,19.2,91.2,123.0,1.56],
 ["254 x 254 x 167",167.1,289.1,265.2,19.2,31.7,12.7,200.3,3.48,10.4,30000.0,9870.0,11.9,6.81,2080.0,744.0,2420.0,1140.0,0.851,8.48,626.0,213.0,1.63],
 ["254 x 254 x 132",132.0,276.3,261.3,15.3,25.3,12.7,200.3,4.36,13.1,22500.0,7530.0,11.6,6.69,1630.0,576.0,1870.0,878.0,0.85,10.3,319.0,168.0,1.19],
 ["254 x 254 x 107",107.1,266.7,258.8,12.8,20.5,12.7,200.3,5.38,15.6,17500.0,5930.0,11.3,6.59,1310.0,458.0,1480.0,697.0,0.848,12.4,172.0,136.0,0.898],
 ["254 x 254 x 89",88.9,260.3,256.3,10.3,17.3,12.7,200.3,6.38,19.4,14300.0,4860.0,11.2,6.55,1100.0,379.0,1220.0,575.0,0.85,14.5,102.0,113.0,0.717],
 ["254 x 254 x 73",73.1,254.1,254.6,8.6,14.2,12.7,200.3,7.77,23.3,11400.0,3910.0,11.1,6.48,898.0,307.0,992.0,465.0,0.849,17.2,57.6,93.1,0.562],
 ["203 x 203 x 127",127.5,241.4,213.9,18.1,30.1,10.2,160.8,2.91,8.88,15400.0,4920.0,9.75,5.5,1280.0,460.0,1520.0,704.0,0.856,7.36,427.0,162.0,0.549],
 ["203 x 203 x 113",113.5,235.0,212.1,16.3,26.9,10.2,160.8,3.26,9.87,13300.0,4290.0,9.59,5.45,1130.0,404.0,1330.0,618.0,0.852,8.11,305.0,145.0,0.464],
 ["203 x 203 x 100",99.6,228.6,210.3,14.5,23.7,10.2,160.8,3.7,11.1,11300.0,3680.0,9.44,5.39,988.0,350.0,1150.0,534.0,0.852,9.01,210.0,127.0,0.386],
 ["203 x 203 x 86",86.1,222.2,209.1,12.7,20.5,10.2,160.8,4.29,12.7,9450.0,3130.0,9.28,5.34,850.0,299.0,977.0,456.0,0.849,10.2,137.0,110.0,0.318],
 ["203 x 203 x 71",71.0,215.8,206.4,10.0,17.3,10.2,160.8,5.09,16.1,7620.0,2540.0,9.18,5.3,706.0,246.0,799.0,374.0,0.853,11.9,80.2,90.4,0.25],
 ["203 x 203 x 60",60.0,209.6,205.8,9.4,14.2,10.2,160.8,6.2,17.1,6120.0,2060.0,8.96,5.2,584.0,201.0,656.0,305.0,0.846,14.1,47.2,76.4,0.197],
 ["203 x 203 x 52",52.0,206.2,204.3,7.9,12.5,10.2,160.8,7.04,20.4,5260.0,1780.0,8.91,5.18,510.0,174.0,567.0,264.0,0.847,15.8,31.8,66.3,0.167],
 ["203 x 203 x 46",46.1,203.2,203.6,7.2,11.0,10.2,160.8,8.0,22.3,4570.0,1550.0,8.82,5.13,450.0,152.0,497.0,231.0,0.847,17.7,22.2,58.7,0.143],
 ["152 x 152 x 51",51.2,170.2,157.4,11.0,15.7,7.6,123.6,4.18,11.2,3230.0,1020.0,7.04,3.96,379.0,130.0,438.0,199.0,0.848,10.1,48.8,65.2,0.061],
 ["152 x 152 x 44",44.0,166.0,155.9,9.5,13.6,7.6,123.6,4.82,13.0,2700.0,860.0,6.94,3.92,326.0,110.0,372.0,169.0,0.847,11.5,31.7,56.1,0.05],
 ["152 x 152 x 37",37.0,161.8,154.4,8.0,11.5,7.6,123.6,5.7,15.5,2210.0,706.0,6.85,3.87,273.0,91.5,309.0,140.0,0.848,13.3,19.2,47.1,0.04],
 ["152 x 152 x 30",30.0,157.6,152.9,6.5,9.4,7.6,123.6,6.98,19.0,1750.0,560.0,6.76,3.83,222.0,73.3,248.0,112.0,0.847,16.1,10.5,38.3,0.031],
 ["152 x 152 x 23",23.0,152.4,152.2,5.8,6.8,7.6,123.6,9.65,21.3,1250.0,400.0,6.54,3.7,164.0,52.6,182.0,80.1,0.842,20.6,4.63,29.2,0.021]
].map(a=>({key:a[0],mass:a[1],D:a[2],B:a[3],tw:a[4],tf:a[5],r:a[6],d:a[7],bT:a[8],dt:a[9],
  Ix:a[10],Iy:a[11],rx:a[12],ry:a[13],Zx:a[14],Zy:a[15],Sx:a[16],Sy:a[17],u:a[18],x:a[19],J:a[20],A:a[21],Iw:a[22]}));

const UCmap=Object.fromEntries(UC.map(s=>[s.key,s]));

const TP385_UC={"356 x 406 x 634":[13700,0.86,38.8,421,34400],"356 x 406 x 551":[9240,0.93,31.1,406,28700],"356 x 406 x 467":[5810,1.04,24.3,390,23300],"356 x 406 x 393":[3550,1.18,18.9,376,18800],"356 x 406 x 340":[2340,1.31,15.5,366,15800],"356 x 406 x 287":[1440,1.49,12.3,356,13000],"356 x 406 x 235":[812,1.75,9.54,346,10300],"356 x 368 x 202":[558,1.82,7.16,326,8240],"356 x 368 x 177":[381,2.03,6.08,321,7110],"356 x 368 x 153":[251,2.3,5.13,316,6060],"356 x 368 x 129":[153,2.66,4.17,312,5020],"305 x 305 x 283":[2030,0.899,6.34,259,9190],"305 x 305 x 240":[1271,1.01,5.03,251,7520],"305 x 305 x 198":[734,1.17,3.88,243,5990],"305 x 305 x 158":[378,1.4,2.87,235,4570],"305 x 305 x 137":[249,1.58,2.39,231,3870],"305 x 305 x 118":[161,1.79,1.98,227,3270],"305 x 305 x 97":[91,2.11,1.56,223,2620],"254 x 254 x 167":[626,0.823,1.63,171,3590],"254 x 254 x 132":[319,0.982,1.19,164,2710],"254 x 254 x 107":[172,1.16,0.899,159,2110],"254 x 254 x 89":[102,1.35,0.717,156,1730],"254 x 254 x 73":[58,1.59,0.563,153,1380],"203 x 203 x 127":[427,0.578,0.549,113,1820],"203 x 203 x 113":[305,0.628,0.464,110,1570],"203 x 203 x 100":[210,0.691,0.386,108,1340],"203 x 203 x 86":[137,0.777,0.318,105,1130],"203 x 203 x 71":[80.2,0.899,0.25,102,914],"203 x 203 x 60":[47.2,1.04,0.197,101,734],"203 x 203 x 52":[31.8,1.17,0.167,98.9,632],"203 x 203 x 46":[22.2,1.29,0.143,97.8,548],"152 x 152 x 51":[48.8,0.568,0.061,60.8,376],"152 x 152 x 44":[31.7,0.639,0.05,59.4,315],"152 x 152 x 37":[19.2,0.734,0.04,58,258],"152 x 152 x 30":[10.5,0.871,0.031,56.6,204],"152 x 152 x 23":[4.63,1.09,0.021,55.4,143],
 // Rows below are not in SCI P385 Appendix A. Derived: IT = J and Iw from the table row above,
 // a = sqrt(E*Iw/(G*IT)) with E = 210000, G = 81000 N/mm2 (m), Wn0 = (D - tf)*B/4 (cm2), Sw1 = (D - tf)*B^2*tf/16 (cm4).
 // Formulas reproduce every tabulated P385 row of this family within 0.5 % (Wn0, Sw1) and 0.9 % (a); see docs/SECTION_LIBRARY.md.
 "356 x 406 x 1299":[94500,0.604,133.1,547,91200],  // derived from table row
 "356 x 406 x 1202":[76300,0.624,114.6,530,81100],  // derived from table row
 "356 x 406 x 1086":[60500,0.642,96.1,504,71500],  // derived from table row
 "356 x 406 x 990":[46900,0.671,81.5,487,62800],  // derived from table row
 "356 x 406 x 900":[36400,0.701,68.9,470,55000],  // derived from table row
 "356 x 406 x 818":[27800,0.739,58.6,456,48300],  // derived from table row
 "356 x 406 x 744":[21400,0.778,50,442,42400],  // derived from table row
 "356 x 406 x 677":[16400,0.824,42.9,430,37500],  // derived from table row
 "356 x 406 x 592":[11400,0.888,34.7,413,31500],  // derived from table row
 "356 x 406 x 509":[7390,0.984,27.6,399,26000],  // derived from table row
 "356 x 406 x 477":[5700,1.04,23.68,397,22400]  // MS, derived
};


/* ==== beam-v03 module 7 ==== */
/* ===========================================================================
   2. SECTION DATA COMMON HELPERS
   Family section tables live in js/sections/*-section-data.js.
   This file normalises the selected family data and provides shared helpers.
   =========================================================================== */
// Normalises the active section (PFC, SHS, or UB) to a single shape so the
// solver/checks/report code can stay section-agnostic wherever the physics
// genuinely is agnostic, and branch only where BS 5950 itself branches
// (classification limits, shear area, LTB applicability, strut curve).
// SCI P385 Appendix A torsional properties (Tables A.1, A.2, A.3, A.7, A.8).
// Open sections: [IT cm4, a m, Iw dm6, Wn0 cm2, Sw1 cm4]; PFC additionally
// [.., Wn2 cm2, Sw1, Sw2, Sw3 cm4, e0 mm, esc mm]; hollow: [IT cm4, Wt cm3].
// Note: P385 IT for open sections includes the Appendix B junction correction,
// so it can differ slightly from the section tables value stored in the main tables.
function tp385For(family,key){
  // P385 A.7 contains 3.2/3.6/6.3 mm walls, not the 3.0/3.5/6.0 mm
  // aliases previously added to this table. Never borrow a thicker wall's Wt.
  if(family==='shs' && (S.shsType==='CF'||[3,3.5,6].includes(Number(key.split('x').at(-1))))) return null;
  if(family==='rhs' && S.rhsType==='CF') return null; // P385 Table A.8 is hot-finished only
  const m=family==='ub'?TP385_UB:family==='uc'?TP385_UC:family==='pfc'?TP385_PFC:family==='shs'?TP385_SHS:family==='rhs'?TP385_RHS:null;
  const r=m&&m[key]; if(!r) return null;
  if(family==='pfc') return {IT:r[0],a:r[1],Iw:r[2],Wn0:r[3],Wn2:r[4],Sw1:r[5],Sw2:r[6],Sw3:r[7],e0:r[8],esc:r[9]};
  if(family==='shs'||family==='rhs') return {IT:r[0],Wt:r[1]};
  return {IT:r[0],a:r[1],Iw:r[2],Wn0:r[3],Sw1:r[4]};
}
function activeSection(){ const s=activeSectionBase(); s.tp=tp385For(S.family,s.key); return s; }
function activeSectionBase(){
  if(S.family==='shs'){
    const map = S.shsType==='CF'? SHS_CFmap : SHS_HFmap;
    const s = map[S.shsKey] || Object.values(map)[0];
    // SCI P363: c = h - 3t for EC3 local buckling, including cold-formed SHS.
    // The legacy BS table stores h - 5t for cold-formed SHS.
    const ratio=S.code==='EC3'? (s.D-3*s.t)/s.t : s.dt;
    return {key:s.key, mass:s.mass, D:s.D, B:s.D, tw:s.t, tf:s.t, r:0, d:ratio*s.t,
      bT:ratio, dt:ratio, Ix:s.I, Iy:s.I, rx:s.r, ry:s.r, Zx:s.Z, Zy:s.Z, Sx:s.S, Sy:s.S,
      u:null, x:null, J:s.J, A:s.A, isBox:true, boxType:S.shsType, kind:'box'};
  }
  if(S.family==='rhs'){
    const cf = S.rhsType==='CF';
    const map = cf? RHS_CFmap : RHSmap;
    const s = map[S.rhsKey] || Object.values(map)[0];
    // SCI P363: c = h - 3t / b - 3t for EC3 local buckling of cold-formed hollow
    // sections (the table stores the BS 5950 h - 5t / b - 5t flat widths), as SHS_CF.
    const bT = (cf && S.code==='EC3')? (s.B-3*s.t)/s.t : s.bT;
    const dt = (cf && S.code==='EC3')? (s.D-3*s.t)/s.t : s.dt;
    return {key:s.key, mass:s.mass, D:s.D, B:s.B, tw:s.t, tf:s.t, r:0, d:dt*s.t,
      bT:bT, dt:dt, Ix:s.Ix, Iy:s.Iy, rx:s.rx, ry:s.ry, Zx:s.Zx, Zy:s.Zy, Sx:s.Sx, Sy:s.Sy,
      u:null, x:null, J:s.J, A:s.A, isBox:true, boxType: cf? 'CF' : 'HF', kind:'box'};
  }
  if(S.family==='ub'){
    const s = UBmap[S.ubKey] || UB[0];
    return Object.assign({isBox:false, boxType:null, kind:'I'}, s);
  }
  if(S.family==='uc'){
    const s = UCmap[S.ucKey] || UC[0];
    return Object.assign({isBox:false, boxType:null, kind:'I'}, s);
  }
  const s = PFCmap[S.sectionKey] || PFC[0];
  return Object.assign({isBox:false, boxType:null, kind:'channel'}, s);
}


function pyFromGrade(grade,t){
  const T=[ [16,40,63,80,100,150],
    {S275:[275,265,255,245,235,225],S355:[355,345,335,325,315,295],S460:[460,440,430,410,400,400]} ];
  const bands=T[0], vals=T[1][grade]||T[1].S275;
  for(let i=0;i<bands.length;i++) if(t<=bands[i]) return vals[i];
  return vals[vals.length-1];
}
const KeByGrade={S275:1.2,S355:1.1,S460:1.0};
function fuFromGrade(grade){ return ({S275:410,S355:470,S460:540})[grade]||410; }
// ---- SCI P385 Appendix C closed-form torsion solution, fork-fork ends ----
// phi, phi', phi'', phi''' for: Case 3 (point torque T at alpha*L), Case 4
// (uniform torque, total T), Case 10 (linear torque 0 -> 2T/L, total T; mirror
// for descending). All hyperbolic terms evaluated in exponential-ratio form to
// avoid catastrophic cancellation at large L/a (P385 Section numerical note).
// Units: L,a in mm; T in N.mm (TOTAL applied torque per P385 p.70 convention);
// GIt in N.mm2. Returns grids over x.
function p385Solve(L,aa,GIt,tqs,nGrid){
  const X=L/aa;
  const e2=u=>Math.exp(-2*Math.max(u,0));
  const den=1-e2(X);
  const ff=(u,v)=>Math.exp(u+v-X)*(1-e2(u))*(1-e2(v))/(2*den);   // sinh(u)sinh(v)/sinh(X)
  const fg=(u,v)=>Math.exp(u+v-X)*(1-e2(u))*(1+e2(v))/(2*den);   // sinh(u)cosh(v)/sinh(X)
  const shR=u=>Math.exp(u-X)*(1-e2(u))/den;                       // sinh(u)/sinh(X)
  const chR=u=>Math.exp(u-X)*(1+e2(u))/den;                       // cosh(u)/sinh(X)
  const H=X/2, dH=1+e2(H);
  const c4=w=>{ const q=Math.abs(H-w); return Math.exp(q-H)*(1+e2(q))/dH; }           // cosh(H-w)/cosh(H)
  const s4=w=>{ const q=Math.abs(H-w); return Math.sign(H-w)*Math.exp(q-H)*(1-e2(q))/dH; } // sinh(H-w)/cosh(H)
  // grid: even spacing + torque application points
  const set=new Set(); const n=nGrid||241;
  for(let i=0;i<=n;i++) set.add(+(L*i/n).toFixed(4));
  tqs.forEach(t=>{ if(t.kind==='point'){ const p=t.alpha*L; [p-0.01,p,p+0.01].forEach(v=>{ if(v>=0&&v<=L) set.add(+v.toFixed(4)); }); } });
  const xs=[...set].sort((p,q)=>p-q);
  const phi=xs.map(()=>0), p1=xs.map(()=>0), p2=xs.map(()=>0), p3=xs.map(()=>0);
  const addCase3=(alpha,T)=>{
    xs.forEach((x,i)=>{
      let al=alpha, xx=x, sgn=1;
      if(x>alpha*L+1e-9){ al=1-alpha; xx=L-x; sgn=-1; }   // mirror for the far segment
      const u1=(1-al)*X, v=xx/aa;
      phi[i]+= (T*aa/GIt)*((1-al)*v - ff(u1,v));
      p1[i] += sgn*(T/GIt)*((1-al) - fg(u1,v));
      p2[i] += -(T/(GIt*aa))*ff(u1,v);
      p3[i] += -sgn*(T/(GIt*aa*aa))*fg(u1,v);
    });
  };
  const addCase4=(T)=>{
    const t=T/L;
    xs.forEach((x,i)=>{
      const w=x/aa;
      phi[i]+= (t*aa*aa/GIt)*( x*(L-x)/(2*aa*aa) + c4(w) - 1 );
      p1[i] += (t*aa/GIt)*( (L-2*x)/(2*aa) - s4(w) );
      p2[i] += (t/GIt)*( c4(w) - 1 );
      p3[i] += -(t/(GIt*aa))*s4(w);
    });
  };
  const addCase10=(T,mirror)=>{
    xs.forEach((x,i)=>{
      const xx=mirror? L-x : x, sgn=mirror? -1 : 1, v=xx/aa;
      phi[i]+= (2*T/GIt)*( xx/6 - xx*aa*aa/(L*L) + (aa*aa/L)*shR(v) - xx*xx*xx/(6*L*L) );
      p1[i] += sgn*(2*T/GIt)*( 1/6 - aa*aa/(L*L) + (aa/L)*chR(v) - xx*xx/(2*L*L) );
      p2[i] += (2*T/GIt)*( shR(v)/L - xx/(L*L) );
      p3[i] += sgn*(2*T/GIt)*( chR(v)/(aa*L) - 1/(L*L) );
    });
  };
  tqs.forEach(t=>{
    if(t.kind==='point') addCase3(Math.min(Math.max(t.alpha,1e-6),1-1e-6),t.T);
    else if(t.kind==='ud') addCase4(t.T);
    else if(t.kind==='lin') addCase10(t.T,!!t.mirror);
  });
  return {xs,phi,p1,p2,p3,X};
}

function ctBoxEN10210(sec){
  // Torsional constants for a hot-finished hollow section per EN 10210-2 Annex A
  // (mean corner radius 1.25t). Verified against SCI P363 section tables It and Ct.
  const t=sec.tf, h=sec.D, b=sec.B, Rc=1.25*t;
  const hp=2*((h-t)+(b-t))-2*Rc*(4-Math.PI);
  const Ah=(h-t)*(b-t)-Rc*Rc*(4-Math.PI);
  const K=2*Ah*t/hp;
  const It=t*t*t*hp/3+2*K*Ah;
  return {Ct:It/(t+K/t), It:It};
}
function defaultRobertson(family,boxType,tf){
  // BS 5950-1:2000 Table 23 (allocation of strut curve) -> Table 24(a..d) Robertson
  // constants a=2.0/3.5/5.5/8.0 for curves a/b/c/d respectively.
  if(family==='shs') return boxType==='CF'? {x:5.5,y:5.5} : {x:2.0,y:2.0}; // HF box: curve a both axes; CF box: curve c both axes
  if(family==='rhs') return boxType==='CF'? {x:5.5,y:5.5} : {x:2.0,y:2.0}; // HF box: curve a both axes; CF box: curve c both axes (Table 23)
  if(family==='ub'){
    // Rolled I-section: x-x/y-y = a/b up to 40mm flange, b/c over 40mm.
    return (tf||0)<=40 ? {x:2.0,y:3.5} : {x:3.5,y:5.5};
  }
  if(family==='uc'){
    // Rolled H-section: x-x/y-y = b/c up to 40mm flange, c/d over 40mm   one
    // curve lower than a rolled I-section at the same thickness (Table 23).
    return (tf||0)<=40 ? {x:3.5,y:5.5} : {x:5.5,y:8.0};
  }
  // PFC (channel): Table 23 directs channels to Table 25 (angle/channel/T-section
  // struts), a distinct method from the generic a/b/c/d curves used here.
  // Channel compression is blocked from PASS unless verified Table 25/section tables
  // data is implemented or supplied.
  return {x:5.5,y:5.5};
}

/* Mirrored PFC: drawing orientation only (web back on the right, toes left).
   All section properties are mirror-invariant; only the shear-centre side,
   the self-weight eccentricity sign, and the drawn geometry flip. */
function pfcMirrored(){ return S.family==='pfc' && !!S.pfcMirror; }
/* Welded bottom plate (continuous, full length, welded to the flange/section
   edges). Modelling + self-weight only - the plate is never designed. */
function plateGeom(sec){
  const p=S.plate;
  if(!p || !p.on) return null;
  const B=+sec.B||0, D=+sec.D||0, t=Math.max(+p.t||0,0);
  const oL=Math.max(+p.outL||0,0), oR=Math.max(+p.outR||0,0);
  if(t<=0) return null;
  const backX=isFinite(+sec.x)? +sec.x : B*0.25;
  const minX=sec.kind==='channel'? (pfcMirrored()? backX-B : -backX) : -B/2;
  const maxX=sec.kind==='channel'? minX+B : B/2;
  const x1=minX-oL, x2=maxX+oR, w=x2-x1;
  const side=p.side==='top'? 'top':'bottom';
  const z1=side==='top'? D/2 : -D/2-t, z2=side==='top'? D/2+t : -D/2;
  return {x1, x2, w, t, oL, oR, side, z1, z2,
    cx:(x1+x2)/2, massPerM:+(w*t*7.85e-3).toFixed(2)};   // 7850 kg/m3 steel
}
function selfWeightValue(sec){
  const pl=plateGeom(sec);
  return +(((+sec.mass||0)+(pl? pl.massPerM:0))*9.81/1000).toFixed(4); // kg/m -> kN/m
}
function selfWeightEccentricity(sec){
  // Self-weight acts through the centroid. For doubly symmetric sections that is
  // also the shear centre; for PFCs the centroid is horizontally eccentric.
  // A one-sided bottom plate shifts the combined centroid - include its share.
  const escBase = sec.kind==='channel' && sec.tp && sec.tp.esc!=null && isFinite(+sec.tp.esc) ? +sec.tp.esc : 0;
  const esc = pfcMirrored()? -escBase : escBase;
  const pl=plateGeom(sec);
  if(!pl || pl.massPerM<=0) return esc;
  const eP=pl.cx+esc;                       // plate centroid offset from the shear centre
  const m=+sec.mass||0, mp=pl.massPerM;
  return +(((m*esc+mp*eP)/(m+mp)).toFixed(2));
}
function syncSelfWeightLoads(){
  // Self-weight is now automatic in comboLoads(). Remove legacy manual
  // self-weight rows so old states or button clicks cannot double-count it.
  S.loads = S.loads.filter(ld=>!ld.isSelfWeight);
}


/* ==== beam-v03 module 8 ==== */
/* ===========================================================================
   3. STATE + INPUT UI
   =========================================================================== */
const CASE_LABELS={G:'Dead (G)',Q:'Imposed (Q)',W:'Wind (W)',E:'Other (E)'};
/* Default combinations (20 Sep 2026 owner decision): ULS 1.35G + 1.5Q (Eq 6.10)
   and SLS 1.0G + 1.0Q, the reference software defaults. The earlier SLS default
   carried Q only (NA 2.23 "variable actions only", G = 0); an owner who wants
   that criterion sets G = 0 under Load Combinations. */
const DEFAULT_COMBOS=[
  {id:'c1', label:'ULS: 1.35G + 1.5Q (Eq 6.10)',  factors:{G:1.35,Q:1.5,W:0,E:0},   sls:false, on:true},
  {id:'s1', label:'SLS: 1.0G + 1.0Q',             factors:{G:1.0, Q:1.0,W:0,E:0},   sls:true,  on:true},
];
/* Section name as shown to the user (pure): the library keys are "200x75x23"
   (PFC), "150x150x6.3" (SHS) and "457 x 191 x 82" (UB / UC / RHS); every
   family is printed "D x B x t" with a spaced " x " (20 Sep 2026 owner
   request: "100  100  6" was confusing). The keys themselves stay the map
   keys and the drop-list option values. */
function sectionDisplayName(key){ return String(key||'').replace(/\s*x\s*/g,' x '); }
/* Span change (20 Sep 2026 owner request): keep the loads, the internal
   hinges and the intermediate LTB restraints inside the member when L changes
   from Lold to Lnew, so that editing the span never leaves "position must be
   within 0 to L m" on screen. Mutates `state` in place and returns it (pure
   otherwise; unit-tested through tests/harness.cjs):
     - a UDL / trapezoidal load whose end x2 sat at Lold (full-span, or
       running to End 2) follows to Lnew, longer or shorter;
     - any x1 / x2 / pos (udl, trap, point, moment) beyond Lnew is clamped to
       Lnew; a UDL / trap whose start is pushed onto its end (the strip lay
       wholly beyond Lnew) is slid inwards keeping its original length where
       the member allows (a 1 m strip stays a 1 m strip; 20 Sep 2026 review:
       restarting at x1 = 0 turned it silently into a full-span load);
     - a restraint position beyond Lnew is clamped to Lnew (kept as the
       string the restraint editor stores; a restraint on the end coincides
       with the end restraint, harmless);
     - a hinge is kept strictly INSIDE the member (validateInputs refuses a
       hinge at a span end): one that would land on or beyond Lnew keeps its
       fraction of the span (7 of 8 m -> 4.375 of 5 m), or, when Lold is
       unknown, moves to 0.9 Lnew (20 Sep 2026 review: the first cut clamped
       it onto the end and the error stayed on screen).
   Loads already inside the new span are untouched and every value stays
   editable afterwards. A non-finite or non-positive Lnew is ignored. */
function clampLoadsToSpan(state, Lold, Lnew){
  state=state||S;
  Lold=+Lold; Lnew=+Lnew;
  if(!(Number.isFinite(Lnew)&&Lnew>0)) return state;
  const oldOk=Number.isFinite(Lold)&&Lold>0;
  const atOldEnd=v=>oldOk&&Math.abs(+v-Lold)<=1e-6;
  const r4=v=>Math.round(v*1e4)/1e4;
  (state.loads||[]).forEach(ld=>{
    if(ld.type==='udl'||ld.type==='trap'){
      const len=Math.max(0,+ld.x2-+ld.x1);   // strip length before the clamp
      if(atOldEnd(ld.x2)||+ld.x2>Lnew) ld.x2=Lnew;
      if(+ld.x1>Lnew) ld.x1=Lnew;
      if(+ld.x1>=+ld.x2 && len>0) ld.x1=r4(Math.max(0,+ld.x2-len));   // slide inwards, same length (a full member at most)
    } else if(ld.type==='point'||ld.type==='moment'){
      if(+ld.pos>Lnew) ld.pos=Lnew;
    }
  });
  (state.hinges||[]).forEach(h=>{
    if(+h.pos>=Lnew-1e-6){
      let p= oldOk? +h.pos*Lnew/Lold : NaN;
      if(!(p>1e-6 && p<Lnew-1e-6)) p=Math.min(+h.pos, 0.9*Lnew);   // Lold unknown, or the hinge sat on / beyond the old end
      h.pos=r4(p);
    }
  });
  (state.ltbRestraints||[]).forEach(r=>{ if(+r.pos>Lnew) r.pos= typeof r.pos==='string'? String(Lnew) : Lnew; });
  return state;
}
/* ---------------------------------------------------------------------------
   Member ends (19 Sep 2026 scope change: the tool handles ONE span, x = 0 =
   End 1 to x = L = End 2). Each end carries six degrees of freedom plus
   warping, each RESTRAINED (true) or FREE (false):
     ux   axial translation                -> axial load path of N_Ed
     uy   lateral translation (minor axis) -> LTB eigen v = 0; minor-axis strut
     uz   vertical translation             -> in-plane solver (vertical support)
     rx   twist about the member axis      -> LTB eigen phi = 0; torsion FE phi = 0
     ry   rotation about the major axis    -> in-plane solver (rotational fixity)
     rz   rotation about the minor axis    -> LTB eigen v' = 0; minor-axis strut
     warp warping (phi' = 0)               -> LTB eigen phi' = 0; torsion FE phi' = 0
   plus, for an end that carries a vertical reaction (uz): holdDown (uplift
   resisted by a designed connection), ss (stiff bearing length, mm; blank =
   the lower bound 0) and stiff (bearing stiffener provided). Internal hinges
   (S.hinges) stay: an in-plane moment release at x, lateral / twist
   continuity kept. Intermediate supports, overhangs and multi-span logic are
   out of scope.
   In-plane mapping (endInPlaneType): uz + ry = 'fixed', uz only = 'pinned',
   ry only = 'guided' (sliding end: rotation held, vertical free), neither =
   'none' (free end, absent from the solver's support list).
   --------------------------------------------------------------------------- */
const END_DOFS=['ux','uy','uz','rx','ry','rz','warp'];
/* Long labels (title / tooltip) and the short labels of the End conditions panel. */
const END_DOF_LABELS={ux:'U<sub>x</sub> axial translation (axial load path of N<sub>Ed</sub>)',uy:'U<sub>y</sub> lateral translation (LTB v = 0; minor-axis strut)',uz:'U<sub>z</sub> vertical translation (vertical support)',rx:'R<sub>x</sub> twist about the member axis (LTB / torsion &phi; = 0)',ry:'R<sub>y</sub> rotation about the major axis (in-plane fixity)',rz:'R<sub>z</sub> rotation about the minor axis (LTB v&prime; = 0: laterally clamped; minor-axis strut)',warp:'warping &phi;&prime; = 0 (LTB eigen and torsion FE)'};
const END_DOF_SHORT={ux:'U<sub>x</sub> axial',uy:'U<sub>y</sub> lateral',uz:'U<sub>z</sub> vertical',rx:'R<sub>x</sub> twist',ry:'R<sub>y</sub> major rotation',rz:'R<sub>z</sub> minor rotation',warp:'Warping &phi;&prime;'};
/* Presets (pure): every preset restrains U_x at End 1 only, so the axial load
   path is statically determinate; tick U_x at End 2 for a fully built-in far
   end (printed as "axial statically indeterminate: N taken as applied").
   Warping is restrained only at a cantilever root (the seventh flag of the
   scope note); every other preset leaves warping free (fork ends).
   `label` is the End conditions drop-list text, `name` the short name printed
   on the End conditions line of the brief and the report. */
const END_PRESETS={
  'ss':           {label:'Simply supported', name:'simply supported', e1:{ux:1,uy:1,uz:1,rx:1}, e2:{uy:1,uz:1,rx:1}},
  'fixed-fixed':  {label:'Fixed-fixed', name:'fixed-fixed', e1:{ux:1,uy:1,uz:1,rx:1,ry:1,rz:1}, e2:{uy:1,uz:1,rx:1,ry:1,rz:1}},
  'fixed-pinned': {label:'Fixed-pinned (propped cantilever)', name:'fixed-pinned (propped cantilever)', e1:{ux:1,uy:1,uz:1,rx:1,ry:1,rz:1}, e2:{uy:1,uz:1,rx:1}},
  'cantilever':   {label:'Cantilever (End 1 fixed, End 2 free)', name:'cantilever', e1:{ux:1,uy:1,uz:1,rx:1,ry:1,rz:1,warp:1}, e2:{}},
  'guided-fixed': {label:'Guided-fixed (End 1 fixed, End 2 sliding: rotation held, vertical free)', name:'guided-fixed', e1:{ux:1,uy:1,uz:1,rx:1,ry:1,rz:1}, e2:{uy:1,rx:1,ry:1,rz:1}},
  'pinned-guided':{label:'Pinned-guided (End 1 pinned, End 2 sliding: rotation held, vertical free)', name:'pinned-guided', e1:{ux:1,uy:1,uz:1,rx:1}, e2:{uy:1,rx:1,ry:1,rz:1}}
};
END_PRESETS['fixed-guided']=END_PRESETS['guided-fixed'];   // alias: End 1 fixed, End 2 guided
/* Normalised flag set of one end (pure; tolerant of a partial object). */
function endFlags(o){
  const e={};
  END_DOFS.forEach(k=>{ e[k]=!!(o&&o[k]); });
  e.holdDown=!!(o&&o.holdDown);
  e.ss=(o&&o.ss!=null&&o.ss!==''&&Number.isFinite(+o.ss))? +o.ss : null;
  e.stiff=!!(o&&o.stiff);
  return e;
}
/* Ends object of a preset with optional per-end overrides {e1:{...}, e2:{...}} (pure). */
function endsPreset(name,overrides){
  const p=END_PRESETS[name];
  if(!p) throw new Error('Unknown end preset "'+name+'" (ss | fixed-fixed | fixed-pinned | cantilever | guided-fixed | pinned-guided)');
  const ov=overrides||{};
  return {e1:endFlags(Object.assign({},p.e1,ov.e1||{})), e2:endFlags(Object.assign({},p.e2,ov.e2||{}))};
}
/* Name of the preset the ends match on the seven DOF flags, else null (pure). */
function endsPresetName(ends){
  const cur=[ends&&ends.e1,ends&&ends.e2].map(endFlags);
  for(const k of Object.keys(END_PRESETS)){
    if(k==='fixed-guided') continue;
    const p=endsPreset(k);
    if([p.e1,p.e2].every((e,i)=>END_DOFS.every(d=>e[d]===cur[i][d]))) return k;
  }
  return null;
}
/* The two ends as a list [{key, n, x (m), ...flags}] (pure; reads S.L). */
function endsList(st){
  st=st||S;
  const E=st.ends||{};
  return [Object.assign({key:'e1',n:1,x:0},endFlags(E.e1)), Object.assign({key:'e2',n:2,x:+st.L},endFlags(E.e2))];
}
function endInPlaneType(e){ return e.uz&&e.ry? 'fixed' : e.uz? 'pinned' : e.ry? 'guided' : 'none'; }
/* In-plane support list the solver and the checks expect: [{end, pos (m),
   type 'pinned'|'fixed'|'guided', holdDown, ss, stiff}]; a free end is absent. */
function endsToSupports(st){
  return endsList(st).map(e=>{ const type=endInPlaneType(e); return type==='none'? null : {end:e.n,pos:e.x,type,holdDown:e.holdDown,ss:e.ss,stiff:e.stiff}; }).filter(Boolean);
}
/* Vertically held ends (a reaction is possible) with their solver-reaction index. */
function verticalEnds(st){ return endsToSupports(st).map((sp,i)=>Object.assign({i},sp)).filter(sp=>sp.type!=='guided'); }
/* In-plane cantilever: End 1 fixed (U_z + R_y), End 2 free of both. */
function isCantilever(st){ const [e1,e2]=endsList(st); return !!(e1.uz&&e1.ry&&!e2.uz&&!e2.ry); }
/* The cantilever NCCI SN006a Tables 3.1-3.3 describe (pure): an in-plane
   cantilever whose root (End 1) also holds U_y, R_z and R_x (v = v' = phi = 0;
   the root warping flag selects the "warping restrained / free" column) and
   whose tip (End 2) is free of every lateral restraint (no U_y, R_x, R_z or
   warping). A cantilever with a laterally / torsionally restrained tip, or a
   root that releases R_x or R_z, is outside the tables (19 Sep 2026 review). */
function isSn006aCantilever(st){
  const [e1,e2]=endsList(st);
  return isCantilever(st) && !!(e1.uy&&e1.rz&&e1.rx) && !(e2.uy||e2.rx||e2.rz||e2.warp);
}
/* An end that is free vertically (cantilever tip, guided tip): the deflection
   limit is L/divisorCant instead of span/divisor. */
function hasFreeVerticalEnd(st){ const [e1,e2]=endsList(st); return !(e1.uz&&e2.uz); }
/* LTB boundary conditions of the two ends for the eigen model: v from U_y,
   v' from R_z, phi from R_x, phi' from warping; an end restraining none of
   them (free end) contributes nothing. Positions in mm. */
function ltbEndRestraints(st){
  return endsList(st).filter(e=>e.uy||e.rz||e.rx||e.warp).map(e=>({x:e.x*1000,v:e.uy?1:0,vp:e.rz?1:0,phi:e.rx?1:0,phip:e.warp?1:0,end:e.n}));
}
/* Warping constant the LTB model uses (dm6): the P385 value where the section
   carries one, else the section tables column, else 0 (hollow sections). Pure. */
function sectionWarpingIw(sec){ const tp=(sec&&sec.tp)||{}; const v=tp.Iw!=null? tp.Iw : (sec? sec.Iw : null); return Number.isFinite(+v)? +v : 0; }
/* Whether a warping flag (phi' = 0) is a boundary condition of the twist
   models for this section. A closed (box) section has I_w = 0: its twist
   equation E I_w phi'''' - G I_T phi'' = m_t reduces to the second-order
   St Venant form, phi' is not a boundary condition, and warping restraint
   may be neglected (EN 1993-1-1 6.2.7(7)). Imposing phi' = 0 on the Hermite
   node of an I_w = 0 model was a spurious constraint whose effect faded with
   the mesh (finding F-E of the 19 Sep 2026 single-span library: RHS-07
   0.51 %, SHS-06 0.42 % "mesh error" at the cantilever root). Pure. */
function warpingApplies(sec){ return sectionWarpingIw(sec)>0; }
/* Printed end boundary conditions of the LTB eigen model (pure; shared by
   the brief (compact: "End 1 v, phi = 0; End 2 free") and the report (with
   the position and the fork / clamped / free label)). With `sec` given, a
   warping flag on an I_w = 0 section is printed as not applied. */
function ltbEndBcText(st, sec, compact){
  const noWarp = sec!=null && !warpingApplies(sec);
  return endsList(st).map(e=>{
    const held=[], warpOn = e.warp && !noWarp;
    if(e.uy) held.push('v'); if(e.rz) held.push('v&prime;'); if(e.rx) held.push('&phi;'); if(warpOn) held.push('&phi;&prime;');
    const note = (e.warp && noWarp) ? (compact? ' [warping not applied: I<sub>w</sub> = 0]' : ' [warping flag not applied: I<sub>w</sub> = 0, St Venant twist only, EN 1993-1-1 6.2.7(7)]') : '';
    if(compact) return 'End '+e.n+' '+(held.length? held.join(', ')+' = 0' : 'free')+note;
    const kind = (e.uy && e.rx && !e.rz && !warpOn) ? ' (fork)' : (!held.length ? ' (free)' : (e.uy && e.rx && e.rz) ? ' (laterally clamped)' : '');
    return 'End '+e.n+' x = '+g(e.x,2)+' m: '+(held.length? held.join(', ')+' = 0' : 'none')+kind+note;
  }).join('; ');
}
/* Twist-restrained ends for the torsion models: phi = 0 from R_x, phi' = 0
   from warping. Positions in mm. */
function twistEnds(st){ return endsList(st).filter(e=>e.rx).map(e=>({end:e.n,pos:e.x*1000,warpFix:!!e.warp})); }
/* EC3 standard-route LTB effective-length factor: the entered L_E/L, else 1.0
   (the closed forms take k = k_w = 1 fork ends over the member length; the
   route refuses ends that are not forks - stdMcrEndsStatus). The BS 5950 path
   takes its L_E from bs5950LtbLength() (Table 13 / 14 by the flags). */
function ltbLeFactor(st){ st=st||S; return (st.leFactor!=null&&st.leFactor!==''&&Number.isFinite(+st.leFactor))? +st.leFactor : 1.0; }
/* Strut effective-length factor of one buckling plane from the end fixities
   (SCI P360 Table 6.2 / BS 5950-1 Table 22 style; pure). t1, r1 = translation
   / rotation held at End 1 in that plane, t2, r2 at End 2. */
function lcrAxisFactor(t1,r1,t2,r2){
  const held=(t1?1:0)+(t2?1:0);
  if(held===2){
    const nr=(r1?1:0)+(r2?1:0);
    if(nr===2) return {K:0.7, text:'held in position at both ends and in direction at both ends: 0.7 L'};
    if(nr===1) return {K:0.85, text:'held in position at both ends, in direction at one end: 0.85 L'};
    return {K:1.0, text:'held in position at both ends, in direction at neither: 1.0 L'};
  }
  if(held===1){
    const heldR= t1? r1 : r2, freeR= t1? r2 : r1;
    if(heldR&&!freeR) return {K:2.0, text:'held in position and direction at one end, the other end free: 2.0 L'};
    if(heldR&&freeR) return {K:1.2, text:'held in position and direction at one end, the other end held in direction but not in position (sway permitted, guided): 1.2 L'};
    if(!heldR&&freeR) return {K:2.0, text:'pinned at the held end, the other end guided (rotation held, sway permitted): 2.0 L (theoretical value; not a P360 Table 6.2 row)'};
    return {K:null, mechanism:true, text:'held in position at one end only and in direction at neither end: sway mechanism'};
  }
  return {K:null, mechanism:true, text:'held in position at neither end: mechanism'};
}
/* Strut effective lengths of both axes from the end fixities (pure): y-y
   (major-axis buckling, in-plane) from U_z / R_y, z-z (minor-axis buckling)
   from U_y / R_z. The entered L_E/L factor overrides both. Returns {Ky, Kz,
   basisY, basisZ, basis, mechanismY, mechanismZ, override, defaultKy,
   defaultKz}; the basis texts carry [verify]. */
function lcrDefaults(st){
  st=st||S;
  const [e1,e2]=endsList(st);
  const y=lcrAxisFactor(e1.uz,e1.ry,e2.uz,e2.ry), z=lcrAxisFactor(e1.uy,e1.rz,e2.uy,e2.rz);
  const override=(st.leFactor!=null&&st.leFactor!==''&&Number.isFinite(+st.leFactor));
  const out={defaultKy:y.K, defaultKz:z.K, mechanismY:!!y.mechanism, mechanismZ:!!z.mechanism, override,
    basisY:'y-y: '+y.text, basisZ:'z-z: '+z.text};
  if(override){
    out.Ky=+st.leFactor; out.Kz=+st.leFactor;
    out.basis='user L<sub>E</sub>/L factor '+(+st.leFactor).toFixed(2)+' on both axes (the end fixities would give '+(y.K!=null? 'y-y '+y.K.toFixed(2) : 'y-y mechanism')+', '+(z.K!=null? 'z-z '+z.K.toFixed(2) : 'z-z mechanism')+') [verify]';
  } else {
    out.Ky=y.K; out.Kz=z.K;
    out.basis='from the end fixities (P360 Table 6.2 / BS 5950 Table 22 style) &mdash; '+out.basisY+'; '+out.basisZ+' [verify]';
  }
  return out;
}
/* ---- End conditions panel state (UI half of the 19 Sep 2026 scope change) ----
   st.endPreset is the drop-list selection: a preset key or 'custom'. The DOF
   flags are the truth; the selection only names them. */
/* Apply a preset to the state: the seven flags of both ends from the preset,
   the seating / hold-down / stiffener entries kept, the hinges cleared (a
   preset can leave an existing hinge as a mechanism). Returns st. */
function applyEndPreset(st,name){
  st=st||S;
  if(name==='fixed-guided') name='guided-fixed';
  const keep=k=>{ const e=(st.ends&&st.ends[k])||{}; return {holdDown:e.holdDown, ss:e.ss, stiff:e.stiff}; };
  st.ends=endsPreset(name,{e1:keep('e1'),e2:keep('e2')});   // throws for an unknown name
  st.endPreset=name;
  st.hinges=[];
  return st;
}
/* Set one DOF flag of one end; the selection becomes 'custom'. Returns st. */
function setEndDof(st,key,dof,on){
  st=st||S;
  if(key!=='e1'&&key!=='e2') throw new Error('setEndDof: end must be e1 or e2, got "'+key+'"');
  if(END_DOFS.indexOf(dof)<0) throw new Error('setEndDof: unknown degree of freedom "'+dof+'" (ux uy uz rx ry rz warp)');
  if(!st.ends) st.ends=endsPreset('ss');
  st.ends[key][dof]=!!on;
  st.endPreset='custom';
  return st;
}
/* Drop-list value for the state (pure): 'custom' when the flags were edited by
   hand, else the preset the flags match, else 'custom'. A stored preset whose
   flags no longer match (a fixture built from endsPreset() alone) is never
   shown: the flags are the truth. */
function endPresetSelection(st){
  st=st||S;
  if(st.endPreset==='custom') return 'custom';
  return endsPresetName(st.ends)||'custom';
}
/* Derived end types (pure): in plane from U_z / R_y, LTB from U_y / R_z / R_x / warping. */
function endInPlaneLabel(t){ return t==='fixed'? 'fixed' : t==='pinned'? 'pinned' : t==='guided'? 'guided (R<sub>y</sub> held, U<sub>z</sub> free)' : 'free'; }
function endLtbType(e){
  if(!(e.uy||e.rx||e.rz||e.warp)) return 'free';
  if(e.uy&&e.rx) return (e.rz? 'laterally clamped' : 'fork')+(e.warp? ', warping fixed' : '');
  const held=[]; if(e.uy) held.push('v'); if(e.rz) held.push('v&prime;'); if(e.rx) held.push('&phi;'); if(e.warp) held.push('&phi;&prime;');
  return 'partial ('+held.join(', ')+' = 0)';
}
/* "End conditions" line of the brief title and the report (pure): the seven
   flags of each end as symbols (W = warping), the derived in-plane and LTB end
   types, the preset name when the flags match one, and the internal hinges. */
function endsConditionsLine(st){
  st=st||S;
  const sym={ux:'U<sub>x</sub>',uy:'U<sub>y</sub>',uz:'U<sub>z</sub>',rx:'R<sub>x</sub>',ry:'R<sub>y</sub>',rz:'R<sub>z</sub>',warp:'W'};
  const parts=endsList(st).map(e=>{
    const held=END_DOFS.filter(k=>e[k]).map(k=>sym[k]);
    return 'End '+e.n+': '+(held.length? held.join(' ')+' restrained' : 'free')+' ('+endInPlaneLabel(endInPlaneType(e))+' in plane; LTB '+endLtbType(e)+')';
  });
  const name=endsPresetName(st.ends);
  const hinges=(st.hinges||[]).length? '; internal hinge'+(st.hinges.length>1? 's' : '')+' at '+st.hinges.map(h=>g(+h.pos,2)+' m').join(', ')+' (in-plane moment release, lateral / twist continuity kept)' : '';
  return parts.join('; ')+' &mdash; '+(name? END_PRESETS[name].name : 'custom end conditions')+hinges;
}
/* Glyph of one end condition for the End conditions panel (pure string
   builder, no DOM): the in-plane symbol (fixed wall, pin or roller, guided
   slider, free end) on a beam stub, the in-plane type above it and the LTB end
   type under it. e = one entry of endsList(); the support sits on the left for
   End 1 and on the right for End 2. */
function endGlyphSvg(e){
  const W=120,H=64, yB=30, sgn=e.n===2? 1 : -1, x0=e.n===2? 88 : 32, xFar=e.n===2? 12 : 108;
  const t=endInPlaneType(e);
  let s='<line x1="'+x0+'" y1="'+yB+'" x2="'+xFar+'" y2="'+yB+'" stroke="#111" stroke-width="3"/>';
  const hatch=(x,y1,y2)=>{ let h=''; for(let y=y1;y<y2;y+=6) h+='<line x1="'+x+'" y1="'+y+'" x2="'+(x+sgn*6)+'" y2="'+(y+6)+'" stroke="#111" stroke-width="1"/>'; return h; };
  if(t==='fixed'){
    s+='<line x1="'+x0+'" y1="'+(yB-16)+'" x2="'+x0+'" y2="'+(yB+16)+'" stroke="#111" stroke-width="2.4"/>'+hatch(x0,yB-16,yB+16);
  } else if(t==='pinned'){
    s+='<polygon points="'+x0+','+yB+' '+(x0-7)+','+(yB+13)+' '+(x0+7)+','+(yB+13)+'" fill="none" stroke="#111" stroke-width="1.6"/>';
    if(e.ux){
      s+='<line x1="'+(x0-11)+'" y1="'+(yB+14)+'" x2="'+(x0+11)+'" y2="'+(yB+14)+'" stroke="#111" stroke-width="1.4"/>';
      for(let x=x0-9;x<=x0+9;x+=6) s+='<line x1="'+x+'" y1="'+(yB+14)+'" x2="'+(x-4)+'" y2="'+(yB+19)+'" stroke="#111" stroke-width="1"/>';
    } else s+='<circle cx="'+(x0-4)+'" cy="'+(yB+16)+'" r="2.2" fill="none" stroke="#111" stroke-width="1.2"/><circle cx="'+(x0+4)+'" cy="'+(yB+16)+'" r="2.2" fill="none" stroke="#111" stroke-width="1.2"/><line x1="'+(x0-11)+'" y1="'+(yB+19)+'" x2="'+(x0+11)+'" y2="'+(yB+19)+'" stroke="#111" stroke-width="1.4"/>';
  } else if(t==='guided'){
    s+='<rect x="'+(x0-3)+'" y="'+(yB-12)+'" width="6" height="24" fill="none" stroke="#111" stroke-width="1.6"/>'+
       '<line x1="'+(x0+sgn*8)+'" y1="'+(yB-16)+'" x2="'+(x0+sgn*8)+'" y2="'+(yB+16)+'" stroke="#111" stroke-width="2.4"/>'+hatch(x0+sgn*8,yB-16,yB+16);
  } else {
    s+='<circle cx="'+x0+'" cy="'+yB+'" r="2.5" fill="#fff" stroke="#111" stroke-width="1.2"/>';
  }
  const top= t==='pinned'? (e.ux? 'pin' : 'roller') : t==='none'? 'free' : t;
  s+='<text x="'+(W/2)+'" y="10" font-family="Arial" font-size="9" fill="#374151" text-anchor="middle">'+top+' in plane</text>';
  const ltb= !(e.uy||e.rx||e.rz||e.warp)? 'free' : (e.uy&&e.rx)? (e.rz? 'clamped' : 'fork')+(e.warp? ' + warping' : '') : 'partial';   // short form of endLtbType() for the glyph width
  s+='<text x="'+(W/2)+'" y="'+(H-4)+'" font-family="Arial" font-size="8.5" fill="#6b7280" text-anchor="middle">LTB: '+ltb+'</text>';
  return '<svg class="end-glyph" viewBox="0 0 '+W+' '+H+'" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="End '+e.n+' condition">'+s+'</svg>';
}



const DEMO={
  code:"EC3",
  family:"ub", sectionKey:"180x75x20", shsType:"HF", shsKey:"150x150x6.3", ubKey:"457 x 191 x 82", ucKey:"203 x 203 x 60", rhsType:"HF", rhsKey:"200 x 100 x 8.0",
  grade:"S275", py:null, anet:null,
  L:8.0,
  ends:endsPreset('ss',{e1:{ss:100},e2:{ss:100}}),   // End 1 (x = 0) / End 2 (x = L) degree-of-freedom flags; ss: stiff bearing length of the seating (mm along the member); blank = lower bound 0
  endPreset:'ss',      // End conditions drop-list: a preset key or 'custom' (the flags above are the truth; endPresetSelection())
  hinges:[],
  loads:[{type:'udl',x1:0,x2:8.0,w:19.7,case:'G'},{type:'udl',x1:0,x2:8.0,w:19.8,case:'Q'}],
  combos:JSON.parse(JSON.stringify(DEFAULT_COMBOS)),
  axial:0, Mz:0, leFactor:null, destab:false, mLTo:null, mxo:null, C1o:null,   // leFactor: blank = strut lengths from the end fixities (lcrDefaults), LTB L_E = L; a number overrides both axes
  LT:null,             // torsional buckling length L_T (m) for a channel under N (cl 6.3.1.4); blank = spacing of the twist restraints (ends with R_x held + restraints with phi held), capped at L_cr,y
  divisor:360, divisorCant:180, deflAbs:null,   // span/360 between two vertically held ends; L/180 when an end is vertically free (cantilever / guided tip, UK NA Table NA.2 [verify]); absolute mm cap (null = none)
  E:210000, Ke:null, robX:null, robY:null,
  restraint:'full',
  mcrMethod:'eigen',   // EC3 unrestrained Mcr: 'eigen' (FE eigensolver, default) | 'standard' (closed form, SN003a/SN006a)
  eccOn:false, za:0,
  ltbRestraints:[],   // intermediate lateral restraints [{pos (m), v, phi, vp, phip}] of the LTB eigen model (js/08-mcr-eigen-patch.js); present here so that a Reset (S = copy of DEMO) keeps the field the patch expects (19 Sep 2026 review)
  zj:0, Cmzo:null,    // eigen patch inputs: monosymmetry z_j (0 for every library section) and the verified C_mz override (null = conservative 1.0; 20 Sep 2026 review: applies to the twist-induced M_z = phi.M_y diagram of EN 1993-6 (A.1) only - with an imposed M_z (constant diagram, psi = 1) C_mz = 1.0 is used whatever the override)
  pfcMirror:false,
  plate:{on:false, side:'bottom', t:10, outL:0, outR:150}
};
let S=JSON.parse(JSON.stringify(DEMO));
function setDesignCode(code){
  const from=S.code;
  if(from===code) return;
  const old=from==='EC3'? [1.35,1.5] : [1.4,1.6];
  const next=code==='EC3'? [1.35,1.5] : [1.4,1.6];
  S.combos.forEach(c=>{
    if(c.id==='c1'&&!c.sls&&c.factors.G===old[0]&&c.factors.Q===old[1]&&c.factors.W===0&&c.factors.E===0){
      c.factors.G=next[0];c.factors.Q=next[1];
      c.label=`ULS: ${next[0]}G + ${next[1]}Q (${code==='EC3'?'Eq 6.10':'BS 5950'})`;
    }
  });
  S.code=code;
  if(S.E===205000||S.E===210000) S.E=code==='EC3'?210000:205000;
}
const $=id=>document.getElementById(id);

/* ---- welded bottom plate UI (modelling + self-weight only) ---- */
function plateState(){ if(!S.plate) S.plate=JSON.parse(JSON.stringify(DEMO.plate)); return S.plate; }
function updatePlateUI(){
  const p=plateState();
  if(!$("platePanel")) return;
  const setV=(id,v)=>{ const el=$(id); if(el && document.activeElement!==el) el.value=v; };
  $("plateOn").checked=!!p.on;
  $("plateFields").style.display=p.on?'':'none';
  if(!p.on) return;
  setV("plateSide",p.side||'bottom');
  setV("plateT",p.t); setV("plateOutL",p.outL); setV("plateOutR",p.outR);
  const sec=activeSection(), pl=plateGeom(sec), h=$("plateHint");
  if(!pl){ h.style.display='none'; return; }
  const sides = pl.oL>0&&pl.oR>0? 'both flange edges' : pl.oL>0? 'the left flange edge' : pl.oR>0? 'the right flange edge' : 'the flange, flush';
  const eSw=selfWeightEccentricity(sec);
  h.innerHTML='<b>Plate '+fmtMM(pl.w)+' &times; '+fmtMM(pl.t)+' mm ('+(pl.side==='top'?'top':'bottom')+')</b> ('+pl.massPerM+' kg/m = '
    +(pl.massPerM*9.81/1000).toFixed(3)+' kN/m), extended past '+sides+'.'
    +' Total self-weight now '+selfWeightValue(sec)+' kN/m'
    +(Math.abs(eSw)>0.5? '; combined self-weight centroid e = '+eSw+' mm from the shear centre.':'.');
  h.style.display='';
}

function fmtMM(v){
  if(!isFinite(v)) return " ";
  const s=(Math.abs(v)<5e-7?0:v).toFixed(1);
  return s.replace(/(\.\d*?)0+$/,'$1').replace(/\.$/,'');
}
function loadHeightReference(sec){
  const depth=+sec.D || 0;
  const wallOrFlange=+sec.tf || 0;
  return {
    topSurface: depth/2,
    topLine: Math.max(depth/2 - wallOrFlange/2, 0),
    bottomSurface: -depth/2
  };
}
function loadHeightReferenceText(sec){
  const ref=loadHeightReference(sec);
  return "Refs from shear centre: top +" + fmtMM(ref.topSurface) + " mm; bottom " +
    fmtMM(ref.bottomSurface) + " mm. Positive z<sub>g</sub> = destabilising.";
}
function loadHeightPerLoadOn(){
  return S.code==='EC3' && (S.restraint||'full')!=='full';
}
function loadZgValue(ld){
  return S.eccOn && ld && ld.zg!=null && isFinite(+ld.zg) ? +ld.zg : (+S.za||0);
}
function syncLoadHeightHint(sec, show){
  const el=$("zaHint");
  if(!el) return;
  el.innerHTML=loadHeightReferenceText(sec);
  el.style.display=show ? '' : 'none';
}

/* HTML of the End conditions panel (pure string builder): two columns, End 1
   (x = 0) and End 2 (x = L), each with its glyph, the seven DOF boxes (checked
   = restrained) and, for an end that holds U_z, the hold-down box and (EC3)
   the stiff bearing length s_s and the bearing-stiffener box. */
function endsPanelHtml(st){
  st=st||S;
  const webBearingOn = st.code==='EC3';
  return endsList(st).map(e=>{
    const dof=k=>`<label class="checkline" title="${END_DOF_LABELS[k].replace(/<[^>]+>/g,'').replace(/"/g,'&quot;')}"><input type="checkbox" data-end="${e.key}" data-dof="${k}"${e[k]?' checked':''}> <span>${END_DOF_SHORT[k]}</span></label>`;
    const extra= e.uz
      ? `<div class="end-extra"><label class="checkline"><input type="checkbox" data-end="${e.key}" data-opt="holdDown"${e.holdDown?' checked':''}> <span>hold-down provided</span></label>${webBearingOn? `
        <div class="fld"><span>Stiff bearing s<sub>s</sub>, mm</span><input type="number" step="1" min="0" placeholder="blank = lower bound 0" value="${e.ss!=null? e.ss : ''}" data-end="${e.key}" data-ss="1" title="Stiff bearing length of the seating along the member (blank = lower bound 0: a station failing at 0 is NOT VERIFIED until the seating length is entered)"></div>
        <label class="checkline"><input type="checkbox" data-end="${e.key}" data-opt="stiff"${e.stiff?' checked':''}> <span>bearing stiffener provided</span></label>` : ''}</div>`
      : `<div class="end-extra end-none">no vertical reaction (U<sub>z</sub> free)</div>`;
    return `<div class="end-col" data-endcol="${e.key}"><div class="end-head">End ${e.n} (x = ${e.n===1? '0' : 'L = '+e.x} m)</div>${endGlyphSvg(e)}
      <div class="end-dofs">${END_DOFS.map(dof).join('')}</div>${extra}</div>`;
  }).join('');
}
function renderEndsPanel(){
  const c=$("endsPanel"); if(!c) return;
  c.innerHTML=endsPanelHtml(S);
  const sel=$("endPreset"); if(sel) sel.value=endPresetSelection(S);
  c.querySelectorAll("[data-dof]").forEach(el=>el.addEventListener("change",e=>{
    setEndDof(S,e.target.dataset.end,e.target.dataset.dof,e.target.checked); renderEndsPanel(); recompute();
  }));
  c.querySelectorAll("[data-opt]").forEach(el=>el.addEventListener("change",e=>{
    S.ends[e.target.dataset.end][e.target.dataset.opt]=e.target.checked; recompute();
  }));
  c.querySelectorAll("[data-ss]").forEach(el=>el.addEventListener("input",e=>{
    // ss: blank = lower bound 0 (NOT VERIFIED if the station fails); a number is kept as entered
    S.ends[e.target.dataset.end].ss = e.target.value===''? null : parseFloat(e.target.value); recompute();
  }));
}
function renderHingeList(){
  const c=$("hingeList"); if(!c) return; c.innerHTML="";
  if(!S.hinges) S.hinges=[];
  S.hinges.forEach((h,i)=>{
    const row=document.createElement("div"); row.className="row";
    row.innerHTML=`<div class="rowhead"><b style="font-size:12px">Internal hinge ${i+1}</b>
      <button class="del" data-hi="${i}">remove</button></div>
      <div class="grid2">
        <div class="fld"><span>Position, m</span><input type="number" step="0.01" value="${h.pos}" data-hp="${i}"></div>
        <div class="fld"><span>Release</span><span style="font-size:11px;color:#6b7280;padding-top:8px">Bending moment M = 0 (in-plane)</span></div>
      </div>`;
    c.appendChild(row);
  });
  c.querySelectorAll("[data-hp]").forEach(el=>el.addEventListener("input",e=>{
    S.hinges[+e.target.dataset.hp].pos=parseFloat(e.target.value); recompute(); }));
  c.querySelectorAll(".del").forEach(b=>b.addEventListener("click",e=>{
    S.hinges.splice(+e.target.dataset.hi,1); renderHingeList(); recompute(); }));
}
function loadFields(ld,i){
  const f=(label,key,val)=>`<div class="fld"><span>${label}</span><input type="number" step="0.01" value="${val}" data-ld="${key}" data-i="${i}"></div>`;
  const caseSel=`<div class="fld"><span>Case</span><select data-ld="case" data-i="${i}">${Object.entries(CASE_LABELS).map(([k,v])=>`<option value="${k}"${ld.case===k?' selected':''}>${v}</option>`).join("")}</select></div>`;
  if(ld.type==='point') return `<div class="grid3">${f("Position, m","pos",ld.pos)}${f("P, kN (?)","P",ld.P)}${caseSel}</div>`;
  if(ld.type==='moment') return `<div class="grid3">${f("Position, m","pos",ld.pos)}${f("M, kN m (?)","M",ld.M)}${caseSel}</div>`;
  if(ld.type==='udl') return `<div class="grid2">${f("Start x1, m","x1",ld.x1)}${f("End x2, m","x2",ld.x2)}</div><div class="grid2" style="margin-top:6px">${f("w, kN/m (?)","w",ld.w)}${caseSel}</div>`;
  if(ld.type==='trap') return `<div class="grid2">${f("Start x1, m","x1",ld.x1)}${f("End x2, m","x2",ld.x2)}</div><div class="grid3" style="margin-top:6px">${f("w1, kN/m","w1",ld.w1)}${f("w2, kN/m","w2",ld.w2)}${caseSel}</div>`;
  return "";
}
/* EN 1993-1-5 clause 6 inputs (EC3 path): per point load and per end a
   stiff bearing length s_s and a "bearing stiffener provided" switch. */
function webBearingInputsOn(){ return S.code==='EC3'; }
function loadBearingFields(ld,i){
  if(!webBearingInputsOn() || ld.type!=='point' || ld.isSelfWeight) return '';
  const ssVal = (ld.ss!=null && ld.ss!=='' && isFinite(+ld.ss))? ld.ss : '';
  return `<div class="grid2" style="margin-top:6px"><div class="fld"><span>Stiff bearing s<sub>s</sub>, mm (blank = 0)</span><input type="number" step="1" min="0" placeholder="0" value="${ssVal}" data-ld="ss" data-i="${i}"></div>
    <label class="checkline" style="align-self:end"><input type="checkbox" data-ldc="stiff" data-i="${i}"${ld.stiff?' checked':''}> <span>bearing stiffener provided (EN 1993-1-5 9.4)</span></label></div>`;
}
function loadOffsetFields(ld,i){
  if(!S.eccOn || ld.type==='moment' || ld.isSelfWeight) return '';
  const eVal = ld.e??0;
  const zgFld = loadHeightPerLoadOn()
    ? `<div class="fld"><span>z<sub>g</sub>, mm (0 = shear centre, + above)</span><input type="number" step="1" value="${loadZgValue(ld)}" data-ld="zg" data-i="${i}"></div>`
    : `<div></div>`;
  return `<div class="grid2" style="margin-top:6px"><div class="fld"><span>e, mm (shear-centre offset; 0 = through shear centre)</span><input type="number" step="1" value="${eVal}" data-ld="e" data-i="${i}"></div>${zgFld}</div>`;
}
function renderLoadList(){
  const c=$("loadList"); c.innerHTML="";
  if(S.loads.length===0) c.innerHTML='<div style="font-size:11px;color:#6b7280;margin:2px 0 6px;">No loads yet &mdash; add one below.</div>';
  S.loads.forEach((ld,i)=>{
    const row=document.createElement("div"); row.className="row";
    const types={point:"Point load",udl:"UDL / partial UDL",trap:"Trapezoidal",moment:"Applied moment"};
    const tag = ld.isSelfWeight? ` <span style="font-size:10px;font-weight:700;color:#7a4;border:1px solid #bcd9a0;background:#f2f8ec;border-radius:4px;padding:1px 5px;">self-weight</span>` : '';
    row.innerHTML=`<div class="rowhead">
      <select data-lt="${i}">${Object.entries(types).map(([k,v])=>`<option value="${k}"${ld.type===k?' selected':''}>${v}</option>`).join("")}</select>${tag}
      <button class="del" data-li="${i}">remove</button></div>${loadFields(ld,i)}${loadOffsetFields(ld,i)}${loadBearingFields(ld,i)}`;
    c.appendChild(row);
  });
  c.querySelectorAll("[data-ld]").forEach(el=>el.addEventListener(el.tagName==='SELECT'?"change":"input",e=>{
    const i=+e.target.dataset.i,k=e.target.dataset.ld;
    // ss: blank = default (0 mm at a point load); a number is kept as entered
    S.loads[i][k]= k==='case'? e.target.value : k==='ss'? (e.target.value===''? null : parseFloat(e.target.value)) : parseFloat(e.target.value);
    recompute();
  }));
  c.querySelectorAll("[data-ldc]").forEach(el=>el.addEventListener("change",e=>{
    const i=+e.target.dataset.i,k=e.target.dataset.ldc;
    S.loads[i][k]=e.target.checked; recompute();
  }));
  c.querySelectorAll("[data-lt]").forEach(sel=>sel.addEventListener("change",e=>{
    const i=+e.target.dataset.lt, t=e.target.value, L=S.L, cs=S.loads[i].case||'Q';
    const old=S.loads[i];
    const base={point:{type:'point',pos:+(L/2).toFixed(3),P:10,case:cs},
                udl:{type:'udl',x1:0,x2:L,w:9,case:cs},
                trap:{type:'trap',x1:0,x2:L,w1:0,w2:12,case:cs},
                moment:{type:'moment',pos:+(L/2).toFixed(3),M:10,case:cs}}[t];
    if(t!=='moment'){ base.e=old.e??0; base.zg=old.zg??(+S.za||0); }
    S.loads[i]=base; renderLoadList(); recompute();
  }));
  c.querySelectorAll(".del").forEach(b=>b.addEventListener("click",e=>{
    S.loads.splice(+e.target.dataset.li,1); renderLoadList(); recompute(); }));
}
function comboLabelFromFactors(combo){
  const parts=['G','Q','W','E'].filter(k=>Math.abs(+combo.factors[k]||0)>1e-9).map(k=>(+combo.factors[k])+k);
  return (combo.sls?'SLS: ':'ULS: ')+(parts.length? parts.join(' + '):'no actions');
}
function renderComboList(){
  const c=$("comboList"); c.innerHTML="";
  S.combos.forEach((combo,i)=>{
    const row=document.createElement("div"); row.className="row";
    const cases=['G','Q','W','E'];
    const factorFlds=cases.map(cs=>`<div class="fld"><span>${cs}</span><input type="number" step="0.05" value="${combo.factors[cs]??0}" data-cf="${cs}" data-ci="${i}"></div>`).join("");
    row.innerHTML=`<div class="rowhead">
      <label style="display:flex;align-items:center;gap:6px;flex:1;font-size:12px;font-weight:700;margin:0">
        <input type="checkbox" style="width:auto" data-con="${i}" ${combo.on?'checked':''}>
        <input type="text" data-clabel="${i}" value="${combo.label}" style="border:none;background:transparent;font-weight:700;padding:2px 0;width:100%">
      </label>
      <span style="font-size:10px;font-weight:700;color:${combo.sls?'#166534':'#1e40af'};border:1px solid ${combo.sls?'#bcd9a0':'#b9c6e6'};background:${combo.sls?'#f2f8ec':'#eef2fb'};border-radius:4px;padding:1px 5px;white-space:nowrap">${combo.sls?'SLS':'ULS'}</span>
      <button class="del" data-cdel="${i}">remove</button></div>
      <div class="grid3" style="grid-template-columns:repeat(4,1fr)">${factorFlds}</div>`;
    c.appendChild(row);
  });
  c.querySelectorAll("[data-con]").forEach(el=>el.addEventListener("change",e=>{
    S.combos[+e.target.dataset.con].on=e.target.checked; recompute(); }));
  c.querySelectorAll("[data-clabel]").forEach(el=>el.addEventListener("input",e=>{
    S.combos[+e.target.dataset.clabel].label=e.target.value; recompute(); }));
  c.querySelectorAll("[data-cf]").forEach(el=>el.addEventListener("input",e=>{
    const i=+e.target.dataset.ci, cs=e.target.dataset.cf;
    S.combos[i].factors[cs]=parseFloat(e.target.value);
    // the combination name always reflects the live factors; editing a factor
    // regenerates it (a hand-typed name lasts until the next factor edit)
    S.combos[i].label=comboLabelFromFactors(S.combos[i]);
    const lab=c.querySelector('[data-clabel="'+i+'"]'); if(lab) lab.value=S.combos[i].label;
    recompute(); }));
  c.querySelectorAll("[data-cdel]").forEach(b=>b.addEventListener("click",e=>{
    S.combos.splice(+e.target.dataset.cdel,1); renderComboList(); recompute(); }));
}
/* Drop-list option text of one library section (pure): "D x B x t FAM (mass
   kg/m)", e.g. "150 x 150 x 6.3 SHS (28.1 kg/m)", "200 x 75 x 23 PFC (23.4
   kg/m)"; the option value stays the raw key (sectionDisplayName). */
function sectionOptionText(s,fam){ return `${sectionDisplayName(s.key)} ${fam} (${s.mass} kg/m)`; }
function syncInputs(){
  const pfcSel=$("pfcSelect"); pfcSel.innerHTML="";
  PFC.forEach(s=>{ const o=document.createElement("option"); o.value=s.key;
    o.textContent=sectionOptionText(s,'PFC'); pfcSel.appendChild(o); });
  pfcSel.value=S.sectionKey;

  const shsSel=$("shsSelect"); shsSel.innerHTML="";
  const shsArr = S.shsType==='CF'? SHS_CF : SHS_HF;
  shsArr.forEach(s=>{ const o=document.createElement("option"); o.value=s.key;
    o.textContent=sectionOptionText(s,'SHS'); shsSel.appendChild(o); });
  if(!(S.shsKey in (S.shsType==='CF'? SHS_CFmap : SHS_HFmap))) S.shsKey = shsArr[0].key;
  shsSel.value=S.shsKey;

  const ubSel=$("ubSelect"); ubSel.innerHTML="";
  UB.forEach(s=>{ const o=document.createElement("option"); o.value=s.key;
    o.textContent=sectionOptionText(s,'UB'); ubSel.appendChild(o); });
  if(!(S.ubKey in UBmap)) S.ubKey = UB[0].key;
  ubSel.value=S.ubKey;

  const ucSel=$("ucSelect"); ucSel.innerHTML="";
  UC.forEach(s=>{ const o=document.createElement("option"); o.value=s.key;
    o.textContent=sectionOptionText(s,'UC'); ucSel.appendChild(o); });
  if(!(S.ucKey in UCmap)) S.ucKey = UC[0].key;
  ucSel.value=S.ucKey;

  const rhsSel=$("rhsSelect"); rhsSel.innerHTML="";
  const rhsArr = S.rhsType==='CF'? RHS_CF : RHS;
  rhsArr.forEach(s=>{ const o=document.createElement("option"); o.value=s.key;
    o.textContent=sectionOptionText(s,'RHS'); rhsSel.appendChild(o); });
  if(!(S.rhsKey in (S.rhsType==='CF'? RHS_CFmap : RHSmap))) S.rhsKey = rhsArr[0].key;
  rhsSel.value=S.rhsKey;

  $("code").value=S.code;
  $("family").value=S.family;
  $("shsType").value=S.shsType;
  if($("rhsType")) $("rhsType").value=S.rhsType||'HF';
  $("pfcRow").style.display = S.family==='pfc'? '' : 'none';
  if($("pfcMirror")) $("pfcMirror").checked = !!S.pfcMirror;
  $("shsRow").style.display = S.family==='shs'? '' : 'none';
  $("shsTypeRow").style.display = S.family==='shs'? '' : 'none';
  if($("rhsTypeRow")) $("rhsTypeRow").style.display = S.family==='rhs'? '' : 'none';
  $("ubRow").style.display = S.family==='ub'? '' : 'none';
  $("ucRow").style.display = S.family==='uc'? '' : 'none';
  $("rhsRow").style.display = S.family==='rhs'? '' : 'none';
  $("restraint").value = S.restraint||'full';
  $("eccOn").checked = !!S.eccOn;
  const sciMode = S.code==='EC3' && (S.restraint||'full')==='full';
  const sec=activeSection();
  $("restraintRow").style.display = S.code==='EC3'? '' : 'none';
  $("bsFactorsRow").style.display = S.code==='EC3'? 'none' : '';
  $("ec3FactorsRow").style.display = (S.code==='EC3' && !sciMode)? '' : 'none';
  $("leRow").style.display = sciMode? 'none' : '';
  $("za").value = S.za||0;
  $("zaRow").style.display = (S.code==='EC3' && !sciMode)? '' : 'none';
  syncLoadHeightHint(sec, S.code==='EC3' && !sciMode);
  if($("mcrMethod")){
    $("mcrMethod").value = S.mcrMethod==='standard'? 'standard' : 'eigen';
    $("mcrMethodRow").style.display = (S.code==='EC3' && !sciMode)? '' : 'none';
  }
  $("destabRow").style.display = sciMode? 'none' : '';
  $("robRow").style.display = S.code==='EC3'? 'none' : '';

  $("grade").value=S.grade;
  $("py").value = S.py!=null? S.py : pyFromGrade(S.grade,sec.tf);
  $("anet").value = S.anet!=null? S.anet : "";
  $("length").value=S.L; $("axial").value=S.axial; if($("Mz")) $("Mz").value=S.Mz;
  $("leFactor").value=S.leFactor??"";
  $("destab").checked=S.destab; $("mLTo").value=S.mLTo??""; $("mxo").value=S.mxo??"";
  $("C1o").value=S.C1o??"";
  if($("LT")){ $("LT").value=S.LT??""; $("LTRow").style.display=(S.code==='EC3' && S.family==='pfc')? '' : 'none'; }
  $("divisor").value=S.divisor; $("E").value=S.E; $("Ke").value=S.Ke??"";
  if($("divisorCant")) $("divisorCant").value=S.divisorCant??180;
  if($("deflAbs")) $("deflAbs").value=S.deflAbs??"";
  const autoRob=defaultRobertson(S.family,sec.boxType,sec.tf);
  $("robertsonX").value = S.robX!=null? S.robX : autoRob.x;
  $("robertsonY").value = S.robY!=null? S.robY : autoRob.y;
  renderEndsPanel(); renderHingeList(); renderLoadList(); renderComboList();
  updatePlateUI();
}


/* ==== beam-v03 module 9 ==== */
/* ===========================================================================
   4. ORCHESTRATION   run analysis + shared checks dispatcher
   =========================================================================== */
/* Effective user loads of ONE combination (pure; reads S.loads). Returns one
   "piece" per load: {i, ld, type, case, e, zg, factor, pos (mm) | x1, x2 (mm)
   with w1, w2 (kN/m, unfactored), P (kN), M (kN.m)}. Every load is present
   (factor 0 when its case is not in the combination) so that the solver's
   x-grid is identical across combinations. The single-span model has no
   pattern combinations (19 Sep 2026 scope change): the companion
   combinations of gammaInfCompanions() are the only generated ones. */
function comboLoadPieces(combo){
  const fac=(combo&&combo.factors)||{};
  const out=[];
  S.loads.forEach((ld,i)=>{
    if(ld.isSelfWeight) return;
    const f=fac[ld.case] ?? 0;
    const base={i,ld,type:ld.type,case:ld.case,e:+(ld.e||0),zg:ld.zg,factor:f};
    if(ld.type==='point'||ld.type==='moment'){ out.push(Object.assign(base,{pos:(+ld.pos)*1000,P:+(ld.P||0),M:+(ld.M||0)})); return; }
    const x1=(+ld.x1)*1000, x2=(+ld.x2)*1000;
    const w1= ld.type==='udl'? +(ld.w||0) : +(ld.w1||0), w2= ld.type==='udl'? +(ld.w||0) : +(ld.w2||0);
    out.push(Object.assign(base,{x1,x2,w1,w2}));
  });
  return out;
}
/* gamma_G,inf companion combinations (19 Sep 2026 review finding): EN 1990
   6.4.3.1(4) - where a permanent action is favourable (a back span holding
   down an overhang, an end span holding down a lifting support) gamma_G,inf
   applies: 1.0 in the STR set B (Table A1.2(B)) and 0.9 in the EQU set A
   (Table A1.2(A)); the single-source rule (6.4.3.1(4) note, A1.3.1(1)) means
   one gamma_G over the whole member, not span by span. For every analysed
   ULS combination whose G factor exceeds 1.0 two
   companions are formed with factors.G replaced by 1.0 and 0.9 (the variable
   factors unchanged). They are solved for their REACTIONS (uplift / hold-down
   design force and the web-bearing reactions); the moment / shear envelopes
   are not extended (printed limitation). Pure. */
function gammaInfCompanions(ulsCombos){
  const out=[];
  ulsCombos.forEach(cb=>{
    const gF=cb.factors.G??0;
    if(!(gF>1+1e-9)) return;
    [[1.0,'B','STR set B','EN 1990 Table A1.2(B)'],[0.9,'A','EQU set A','EN 1990 Table A1.2(A)']].forEach(([gInf,set,setName,ref])=>{
      out.push(Object.assign({},cb,{id:(cb.id||'combo')+'#gInf'+set,label:cb.label+' [&gamma;<sub>G,inf</sub> = '+gInf.toFixed(1)+', '+setName+']',
        factors:Object.assign({},cb.factors,{G:gInf}),parent:cb,gInf,gInfSet:set,gInfRef:ref,companion:true}));
    });
  });
  return out;
}

function comboLoads(combo){
  // Always include every load (factor 0 if its case isn't in this combo) so that
  // load/support positions   and therefore the solver's x-grid   are identical
  // across every combination. That's what makes the envelope comparison below valid.
  const loads = comboLoadPieces(combo).map(p=>{
    if(p.type==='point') return {type:'point',pos:p.pos,P:-p.P*p.factor*1000};
    if(p.type==='moment') return {type:'moment',pos:p.pos,M:p.M*p.factor*1e6};
    return {type:'udl',x1:p.x1,x2:p.x2,w1:-p.w1*p.factor,w2:-p.w2*p.factor};
  });
  const sw=selfWeightValue(activeSection()), gFactor=combo.factors.G ?? 0;
  if(gFactor!==0){
    loads.push({type:'udl',x1:0,x2:S.L*1000,w1:-sw*gFactor,w2:-sw*gFactor,isAutoSelfWeight:true});
  }
  return loads;
}
function comboHasServiceLoad(combo){
  const eps=1e-12;
  for(const ld of S.loads.filter(ld=>!ld.isSelfWeight)){
    const factor=Math.abs(combo.factors[ld.case] ?? 0);
    if(factor<eps) continue;
    if(ld.type==='point' && Math.abs(ld.P||0)>eps) return true;
    if(ld.type==='moment' && Math.abs(ld.M||0)>eps) return true;
    if(ld.type==='udl' && Math.abs(ld.w||0)>eps) return true;
    if(ld.type==='trap' && (Math.abs(ld.w1||0)>eps || Math.abs(ld.w2||0)>eps)) return true;
  }
  return Math.abs(combo.factors.G ?? 0)>eps && selfWeightValue(activeSection())>eps;
}
/* ---------------------------------------------------------------------------
   Stability of the end-restraint model (pure; reads S.ends, S.hinges, S.L,
   S.axial, S.eccOn, the LTB mode). Every unknown combination throws here with
   a clear message - nothing downstream may return NaN.
   In-plane: the solver's own restraint-rank test (beamRestraintRank) on the
   U_z / R_y flags plus the hinges. Out-of-plane (LTB eigen): at least one end
   with U_y and R_x held; when only one end holds U_y the member is a lateral
   cantilever, allowed only when that end also holds R_z and R_x (its warping
   condition is stated in the report). Torsion: any torque needs R_x at one
   end at least. Axial: N_Ed != 0 needs U_x at one end; both ends held is
   allowed and printed as statically indeterminate. Strut lengths: an axis
   whose fixities form a sway mechanism cannot carry N_Ed.
   Returns {errors:[...], notes:[...]} (notes are printed, not blocking).
   --------------------------------------------------------------------------- */
function endsStability(st,sec){
  st=st||S;
  const errs=[], notes=[];
  const [e1,e2]=endsList(st);
  const L=+st.L;
  const hinges=(st.hinges||[]).map(h=>(+h.pos)*1000).filter(x=>Number.isFinite(x)&&x>1e-6&&x<L*1000-1e-6);
  const sup=endsToSupports(st).map(sp=>({pos:sp.pos*1000,type:sp.type}));
  // in-plane
  if(!e1.uz && !e2.uz) errs.push('Neither end restrains vertical translation U<sub>z</sub>: mechanism (rigid-body vertical translation). Restrain U<sub>z</sub> at one end at least.');
  else if(!(e1.uz&&e2.uz) && !e1.ry && !e2.ry){
    const held=e1.uz? 1 : 2, free=e1.uz? 2 : 1;
    errs.push('End '+free+' has no vertical restraint and End '+held+' does not restrain rotation: mechanism (rigid-body rotation about End '+held+'). Restrain R<sub>y</sub> at End '+held+' (cantilever / propped layout) or U<sub>z</sub> at End '+free+'.');
  }
  if(!errs.length && !beamRestraintRank(L*1000,sup,hinges).ok){
    errs.push(hinges.length
      ? 'Under-restrained layout (mechanism): the '+hinges.length+' internal hinge(s) release the in-plane moment and the end restraints U<sub>z</sub> / R<sub>y</sub> do not hold every segment they separate (each hinge needs one more restraint unit: '+(2+hinges.length)+' are needed, U<sub>z</sub> and R<sub>y</sub> counting one each; e.g. one hinge in a fixed - pinned, fixed - guided or fixed - fixed member, two hinges in a fixed - fixed member). Restrain R<sub>y</sub> at an end or remove a hinge.'
      : 'Under-restrained layout (mechanism): the end restraints U<sub>z</sub> / R<sub>y</sub> do not hold the member in its plane.');
  }
  if(e2.uz&&e2.ry&&!e1.uz&&!e1.ry) errs.push('A cantilever must have its root at End 1 (x = 0) and its free tip at End 2 (x = L): mirror the member.');
  // out-of-plane: whenever an LTB check is evaluated - the EC3 unrestrained
  // routes (eigen / standard) and the BS 5950 path, which always checks LTB
  // (19 Sep 2026 review: the BS path let a released U_y / R_x pass unvalidated)
  const ec3=st.code==='EC3';
  const ltbOn = !ec3 || (st.restraint||'full')!=='full';
  if(ltbOn){
    const forks=[e1,e2].filter(e=>e.uy&&e.rx);
    const lateral=[e1,e2].filter(e=>e.uy);
    const singular= ec3? 'the lateral stiffness matrix of the LTB model is singular' : 'the member has no lateral support for the LTB check (lateral mechanism)';
    if(!lateral.length) errs.push('Lateral-torsional buckling: neither end restrains lateral translation U<sub>y</sub>; '+singular+'. Restrain U<sub>y</sub> (and R<sub>x</sub>) at one end at least'+(ec3? ', or set the member fully restrained.' : '.'));
    else if(!forks.length) errs.push('Lateral-torsional buckling: no end restrains both U<sub>y</sub> and R<sub>x</sub> (a fork / torsional restraint); the twist mode is unrestrained. Restrain R<sub>x</sub> at an end that holds U<sub>y</sub>.');
    else if(lateral.length===1){
      const e=lateral[0];
      if(!(e.rz&&e.rx)) errs.push('Lateral-torsional buckling: only End '+e.n+' restrains lateral translation U<sub>y</sub>, so the member is a lateral cantilever; that end must also restrain R<sub>z</sub> (lateral bending) and R<sub>x</sub> (twist), or '+singular+'.');
      else notes.push('Lateral cantilever: End '+e.n+' is the only end holding U<sub>y</sub>; it restrains R<sub>z</sub> and R<sub>x</sub> and its warping is '+(e.warp? (sec&&!warpingApplies(sec)? 'flagged restrained but not applied (I<sub>w</sub> = 0: &phi;&prime; is not a boundary condition of a closed section, EN 1993-1-1 6.2.7(7))' : 'restrained (&phi;&prime; = 0)') : 'free')+'; the other end is laterally free.');
    }
    if(!(e1.rx||e2.rx)) errs.push('Lateral-torsional buckling: neither end restrains twist R<sub>x</sub>: torsional mechanism.');
  }
  // torsion: any torque needs a twist restraint
  const anyTorque = !!st.eccOn && (st.loads.some(ld=>!ld.isSelfWeight && ld.type!=='moment' && Math.abs(ld.e||0)>1e-9) || (typeof selfWeightEccentricity==='function' && Math.abs(selfWeightEccentricity(activeSection()))>1e-9));
  if(anyTorque && !(e1.rx||e2.rx)) errs.push('Torsion: the loads apply a torque but neither end restrains twist R<sub>x</sub>: torsional mechanism. Restrain R<sub>x</sub> at one end at least.');
  // axial
  const N=+st.axial||0;
  if(Math.abs(N)>1e-9){
    if(!(e1.ux||e2.ux)) errs.push('Axial force N<sub>Ed</sub> = '+N+' kN is entered but neither end restrains axial translation U<sub>x</sub>: no axial load path (mechanism). Restrain U<sub>x</sub> at one end.');
    else if(e1.ux&&e2.ux) notes.push('Axial: both ends restrain U<sub>x</sub> - axial statically indeterminate: N taken as applied (N<sub>Ed</sub> = '+N+' kN over the whole member).');
    const lcr=lcrDefaults(st);
    if(N>0 && !lcr.override){
      if(lcr.mechanismY) errs.push('Strut buckling y-y: the U<sub>z</sub> / R<sub>y</sub> fixities form a mechanism ('+lcr.basisY+'); N<sub>Ed</sub> cannot be carried.');
      if(lcr.mechanismZ) errs.push('Strut buckling z-z: the U<sub>y</sub> / R<sub>z</sub> fixities form a sway mechanism ('+lcr.basisZ+'); restrain U<sub>y</sub> at the other end or R<sub>z</sub> at the held end, or enter an L<sub>E</sub>/L factor.');
    }
  }
  return {errors:errs, notes};
}
function validateInputs(py,E,ulsCombos,slsCombos){
  const errs=[];
  const finite=(v)=>v!==null && v!=='' && Number.isFinite(+v);
  const inSpan=(x)=>finite(x) && +x>=-1e-9 && +x<=S.L+1e-9;
  if(!(finite(S.L) && S.L>0)) errs.push("Member length L must be greater than 0.");
  if(!(finite(py) && py>0)) errs.push("Design strength must be greater than 0.");
  if(!(finite(E) && E>0)) errs.push("E must be greater than 0.");
  if(!(finite(S.divisor) && S.divisor>0)) errs.push("Deflection divisor must be greater than 0.");
  if(S.divisorCant!=null && !(finite(S.divisorCant) && S.divisorCant>0)) errs.push("Cantilever deflection divisor must be greater than 0.");
  if(S.deflAbs!=null && S.deflAbs!=='' && !(finite(S.deflAbs) && S.deflAbs>0)) errs.push("Absolute deflection limit must be blank or greater than 0 mm.");
  if(S.leFactor!=null && S.leFactor!=='' && !(finite(S.leFactor) && S.leFactor>0)) errs.push("Effective length factor must be blank (from the end fixities) or greater than 0.");
  ['axial','Mz','za'].forEach(k=>{ if(!finite(S[k])) errs.push(`${k} must be a finite number.`); });
  const area=activeSection().A;
  if(S.anet!=null && !(finite(S.anet)&&S.anet>0&&S.anet<=area)) errs.push('Net area must be greater than zero and no greater than the gross area.');
  ['Ke','robX','robY','C1o'].forEach(k=>{ if(S[k]!=null && !(finite(S[k])&&S[k]>0)) errs.push(`${k} override must be greater than zero.`); });
  if(S.LT!=null && S.LT!=='' && !(finite(S.LT) && +S.LT>0)) errs.push('Torsional buckling length L_T must be blank (= spacing of the twist restraints) or greater than 0 m.');
  if(S.mcrMethod!=null && S.mcrMethod!=='eigen' && S.mcrMethod!=='standard') errs.push(`Unknown Mcr method "${S.mcrMethod}": use "eigen" (FE eigensolver) or "standard" (closed form).`);
  [['mLTo',0.44],['mxo',0.4]].forEach(([k,min])=>{ if(S[k]!=null && !(finite(S[k])&&S[k]>=min&&S[k]<=1)) errs.push(`${k} override must be between ${min} and 1.`); });
  if(S.mLTo!=null && (S.destab || isCantilever(S)) && S.mLTo!==1) errs.push('mLT must be 1 for cantilevers and destabilising loading.');
  if(!S.ends || !S.ends.e1 || !S.ends.e2) errs.push('The member needs its two ends (S.ends.e1 / S.ends.e2): use a preset (endsPreset) or set the degree-of-freedom flags.');
  else {
    endsList().forEach(e=>{
      END_DOFS.forEach(k=>{ const v=(S.ends[e.key]||{})[k]; if(v!=null && typeof v!=='boolean' && v!==0 && v!==1) errs.push(`End ${e.n} ${k} must be true (restrained) or false (free).`); });
      if(e.ss!=null && !(finite(e.ss) && +e.ss>=0)) errs.push(`End ${e.n} stiff bearing length s_s must be blank (default) or a number >= 0 mm.`);
      const raw=(S.ends[e.key]||{}).ss; if(raw!=null && raw!=='' && !Number.isFinite(+raw)) errs.push(`End ${e.n} stiff bearing length s_s must be blank (default) or a number >= 0 mm.`);
    });
  }
  (S.hinges||[]).forEach((h,i)=>{
    if(!inSpan(h.pos)) errs.push(`Internal hinge ${i+1} position must be within 0 to ${S.L} m.`);
    else if(+h.pos<=1e-6 || +h.pos>=S.L-1e-6) errs.push(`Internal hinge ${i+1} must be inside the span, not at an end (release the end rotation R_y instead).`);
  });
  S.loads.forEach((ld,i)=>{
    const tag=`Load ${i+1}`;
    ['e','zg'].forEach(k=>{ if(ld[k]!=null&&!finite(ld[k])) errs.push(`${tag} ${k} must be a finite number.`); });
    if(ld.ss!=null && ld.ss!=='' && !(finite(ld.ss) && +ld.ss>=0)) errs.push(`${tag} stiff bearing length s_s must be blank (default 0) or a number >= 0 mm.`);
    if(!CASE_LABELS[ld.case]) errs.push(`${tag} has an unknown load case.`);
    if(ld.type==='point'){
      if(!inSpan(ld.pos)) errs.push(`${tag} point-load position must be within 0 to ${S.L} m.`);
      if(!finite(ld.P)) errs.push(`${tag} point load must be numeric.`);
    } else if(ld.type==='moment'){
      if(!inSpan(ld.pos)) errs.push(`${tag} moment position must be within 0 to ${S.L} m.`);
      if(!finite(ld.M)) errs.push(`${tag} moment must be numeric.`);
    } else if(ld.type==='udl' || ld.type==='trap'){
      if(!inSpan(ld.x1) || !inSpan(ld.x2)) errs.push(`${tag} load extents must be within 0 to ${S.L} m.`);
      if(!(finite(ld.x1) && finite(ld.x2) && +ld.x2>+ld.x1)) errs.push(`${tag} must have x2 greater than x1.`);
      if(ld.type==='udl' && !finite(ld.w)) errs.push(`${tag} UDL intensity must be numeric.`);
      if(ld.type==='trap' && (!finite(ld.w1) || !finite(ld.w2))) errs.push(`${tag} trapezoidal intensities must be numeric.`);
    } else errs.push(`${tag} has an unknown load type.`);
  });
  [...ulsCombos,...slsCombos].forEach(combo=>{
    ['G','Q','W','E'].forEach(cs=>{
      if(!finite(combo.factors[cs])) errs.push(`Combination "${combo.label}" has a non-numeric ${cs} factor.`);
    });
  });
  S.loads.filter(ld=>!ld.isSelfWeight).forEach((ld,i)=>{
    const active=ld.type==='point'?Math.abs(ld.P)>0:ld.type==='moment'?Math.abs(ld.M)>0:ld.type==='udl'?Math.abs(ld.w)>0:Math.abs(ld.w1)>0||Math.abs(ld.w2)>0;
    if(active&&!ulsCombos.some(cb=>Math.abs(cb.factors[ld.case])>1e-12)) errs.push(`Load ${i+1} (${ld.case}) is omitted from every enabled ULS combination.`);
  });
  if(!slsCombos.some(comboHasServiceLoad)){
    errs.push("No SLS loads applied: every enabled SLS combination has zero factors for the active load cases. Enable a non-zero SLS factor for a load case that is present, or add a serviceability load.");
  }
  if(!errs.length && S.ends && S.ends.e1 && S.ends.e2) errs.push(...endsStability(S).errors);
  if(errs.length) throw errs.join(" ");
}
/* Memo of warping-torsion FE solves (19 Sep 2026 review, performance): keyed
   by the full input of warpingTorsionFE (L, E I_w, G I_T, support list with
   the warping flags, torque list), bounded; a re-render after an unrelated
   input, or a pattern / companion whose torque list repeats, is served from
   the cache. Pure function results are immutable here (never mutated). */
const TORSION_FE_CACHE=new Map(), TORSION_FE_CACHE_MAX=64;
function warpingTorsionFEMemo(opts,stats){
  const key=JSON.stringify([opts.L,opts.EIw,opts.GIt,opts.supports,opts.torques]);
  const hit=TORSION_FE_CACHE.get(key);
  if(hit){ TORSION_FE_CACHE.delete(key); TORSION_FE_CACHE.set(key,hit); if(stats) stats.cached++; return hit; }
  const sol=warpingTorsionFE(opts);
  if(stats) stats.solved++;
  TORSION_FE_CACHE.set(key,sol);
  if(TORSION_FE_CACHE.size>TORSION_FE_CACHE_MAX) TORSION_FE_CACHE.delete(TORSION_FE_CACHE.keys().next().value);
  return sol;
}
function analyse(){
  const sec=activeSection();
  syncSelfWeightLoads();
  const L=S.L*1000;
  const py = S.py!=null? S.py : pyFromGrade(S.grade,sec.tf);
  const E=S.E;
  const Ix=sec.Ix*1e4, EI=E*Ix;
  const ulsUser=S.combos.filter(c=>c.on && !c.sls);
  const slsUser=S.combos.filter(c=>c.on && c.sls);
  if(ulsUser.length===0) throw 'Enable at least one ULS load combination (see "Load Combinations").';
  if(slsUser.length===0) throw 'Enable at least one SLS (deflection) load combination (see "Load Combinations").';
  validateInputs(py,E,ulsUser,slsUser);
  // In-plane support list of the two ends (compatibility shim of the ends model:
  // uz + ry = fixed, uz = pinned, ry = guided, neither = free / absent)
  const supportsMM=endsToSupports(S).map(s=>({pos:(+s.pos)*1000,type:s.type,end:s.end}));
  const hingesMM=(S.hinges||[]).map(h=>(+h.pos)*1000).filter(x=>x>1e-6 && x<L-1e-6);
  const stability=endsStability(S,sec);
  const ends=endsList(S);
  const cant=isCantilever(S);
  const ulsCombos=ulsUser.slice();
  const slsCombos=slsUser.slice();
  const companionCombos=gammaInfCompanions(ulsCombos);
  const companionNote= companionCombos.length
    ? '&gamma;<sub>G,inf</sub> companions: '+companionCombos.length+' ULS combination(s) re-solved with G at 1.0 (STR set B, EN 1990 6.4.3.1(4) / Table A1.2(B)) and at 0.9 (EQU set A, Table A1.2(A)) over the whole member (single-source permanent action) for the END REACTIONS only - uplift / hold-down forces and the web-bearing reactions; the moment and shear envelopes keep the entered &gamma;<sub>G</sub>, so add a reduced-G combination by hand where a relieving permanent action could increase a moment.'
    : '&gamma;<sub>G,inf</sub> companions: none generated (no enabled ULS combination has a G factor above 1.0; a relieving-G case must be entered by hand where a permanent action is favourable).';

  // Run every enabled ULS combination; the same load/support geometry means every
  // combo's result lands on an identical x-grid, so elementwise envelopes are valid.
  const mechanismMsg=()=> hingesMM.length? "Under-restrained layout (mechanism): an internal hinge has left part of the member unrestrained. Restrain R_y at both ends or remove the hinge." : "Under-restrained layout (mechanism): the end restraints U_z / R_y do not hold the member in its plane.";
  const solveUls=combo=>{
    const loads=comboLoads(combo);
    const r=solveBeam(L,EI,supportsMM,loads,120,hingesMM);
    if(!r.w.every(Number.isFinite)) throw mechanismMsg();
    const fb=sfdBmd(L,supportsMM,loads,r.reactions);
    let Vmax=0; fb.V.forEach(v=>{ if(Math.abs(v)>Math.abs(Vmax)) Vmax=v; });
    let Mmax=0,Mpos=0; fb.xs.forEach((x,i)=>{ if(Math.abs(fb.M[i])>Math.abs(Mmax)){Mmax=fb.M[i];Mpos=x;} });
    return {combo,r,fb,Vmax,Mmax,Mpos};
  };
  const ulsResults=ulsCombos.map(solveUls);
  // gamma_G,inf companions (1.0 STR set B and 0.9 EQU set A of every ULS
  // combination with G > 1.0): reactions for uplift / hold-down and web bearing
  const ulsCompanions=companionCombos.map(solveUls);
  const xs=ulsResults[0].fb.xs;
  const Venv=xs.map((_,i)=>{ let best=0; ulsResults.forEach(res=>{ if(Math.abs(res.fb.V[i])>Math.abs(best)) best=res.fb.V[i]; }); return best; });
  const Menv=xs.map((_,i)=>{ let best=0; ulsResults.forEach(res=>{ if(Math.abs(res.fb.M[i])>Math.abs(best)) best=res.fb.M[i]; }); return best; });
  let governV=ulsResults[0]; ulsResults.forEach(r=>{ if(Math.abs(r.Vmax)>Math.abs(governV.Vmax)) governV=r; });
  let governM=ulsResults[0]; ulsResults.forEach(r=>{ if(Math.abs(r.Mmax)>Math.abs(governM.Mmax)) governM=r; });
  const Vmax=governV.Vmax, Mmax=governM.Mmax, Mpos=governM.Mpos;
  // m-factor inputs (quarter-point moments) must come from ONE moment-diagram shape  
  // the governing-moment combo's own BMD   not a mix of different combos' diagrams.
  const gfb=governM.fb;
  const Mq=interpAt(gfb.xs,gfb.M,L*0.25), Mh=interpAt(gfb.xs,gfb.M,L*0.5), Mq3=interpAt(gfb.xs,gfb.M,L*0.75);
  let M24=0; gfb.xs.forEach((x,i)=>{ if(x>=L*0.25-1&&x<=L*0.75+1) M24=Math.max(M24,Math.abs(gfb.M[i])); });
  // end moments read a fraction inside the member (as analysisForCombination does):
  // the grid closes to zero at x = L beyond an end couple, and a fixed-end sample
  // exactly at the support is the far side of the reaction jump
  const M0end=interpAt(gfb.xs,gfb.M,1e-4), MLend=interpAt(gfb.xs,gfb.M,L-1e-4);
  const reactions=governM.r.reactions;

  // SLS deflection: worst of every enabled SLS combination, over the one span
  // against its limit: span/S.divisor when both ends are held vertically,
  // L/S.divisorCant when an end is vertically free (cantilever tip, guided
  // tip; UK NA to EN 1993-1-1 Table NA.2 cantilever row, default 180
  // [verify]), capped by the optional absolute limit S.deflAbs (mm). A held
  // end has w = 0, so a tip value is the deflection relative to the root.
  const deflCant=hasFreeVerticalEnd(S);
  const divisorCant=(S.divisorCant!=null && Number.isFinite(+S.divisorCant) && +S.divisorCant>0)? +S.divisorCant : 180;
  const deflAbs=(S.deflAbs!=null && S.deflAbs!=='' && Number.isFinite(+S.deflAbs) && +S.deflAbs>0)? +S.deflAbs : null;
  const slsResults=slsCombos.map(combo=>{
    const loads=comboLoads(combo);
    const r=solveBeam(L,EI,supportsMM,loads,120,hingesMM);
    if(!r.w.every(Number.isFinite)) throw mechanismMsg();
    let dmax=0,dpos=0; r.nodes.forEach((x,i)=>{ if(Math.abs(r.w[i])>Math.abs(dmax)){dmax=r.w[i];dpos=x;} });
    const divisor= deflCant? divisorCant : S.divisor;
    const limSpan=L/divisor;
    const absGoverns= deflAbs!=null && deflAbs<limSpan;
    const limit= absGoverns? deflAbs : limSpan;
    const util=Math.abs(dmax)/limit;
    const deflection={no:1,start:0,end:L,span:L,cant:deflCant,dmax,dpos,limit,util,divisor,limSpan,abs:deflAbs,absGoverns,combo:combo.label};
    return {combo,r,dmax,dpos,deflection,segs:[deflection]};
  });
  let governD=slsResults[0]; slsResults.forEach(r=>{ if(r.deflection.util>governD.deflection.util) governD=r; });
  const dmax=governD.dmax, dpos=governD.dpos;
  const deflSegments=[governD.deflection];

  // ---- uplift / hold-down (EN 1990 2.4.4 EQU): every combination's reactions ----
  // A negative vertical reaction (up = positive) means the end must hold the
  // member down. Recorded per combination for every vertically held end; the
  // worst per end is kept (design force = the worst ULS value incl. the
  // gamma_G,inf companions, SLS uplift listed separately). n = end number.
  const vEnds=verticalEnds(S);
  const uplift=[];
  ulsResults.concat(ulsCompanions).forEach(res=>vEnds.forEach(ve=>{ const re=res.r.reactions[ve.i]; if(re.V< -1) uplift.push({n:ve.end,pos:re.pos,R:re.V/1000,combo:res.combo.label,sls:false,gInf:res.combo.gInf||null}); }));
  slsResults.forEach(res=>vEnds.forEach(ve=>{ const re=res.r.reactions[ve.i]; if(re.V< -1) uplift.push({n:ve.end,pos:re.pos,R:re.V/1000,combo:res.combo.label,sls:true}); }));
  const upliftSupports=vEnds.map(ve=>{
    const rows=uplift.filter(u=>u.n===ve.end);
    if(!rows.length) return null;
    const worst=rows.reduce((p,u)=>u.R<p.R?u:p);
    const ulsRows=rows.filter(u=>!u.sls), slsRows=rows.filter(u=>u.sls);
    const worstUls=ulsRows.length? ulsRows.reduce((p,u)=>u.R<p.R?u:p) : null;
    const worstSls=slsRows.length? slsRows.reduce((p,u)=>u.R<p.R?u:p) : null;
    return {n:ve.end,pos:worst.pos,type:ve.type,holdDown:!!ve.holdDown,R:worst.R,combo:worst.combo,sls:worst.sls,
      RUls:worstUls? worstUls.R : null, comboUls:worstUls? worstUls.combo : null, gInfUls:worstUls? worstUls.gInf : null,
      RSls:worstSls? worstSls.R : null, comboSls:worstSls? worstSls.combo : null, nCombos:rows.length};
  }).filter(Boolean);

  // ---- torsion from load eccentricity (loads at e from the shear centre) ----
  // Torque loads mirror the transverse loads: q_T(x) = w(x)*e, point torques P*e.
  // Twist is prevented at every end whose R_x is restrained: GIt*phi'' = -q_T
  // with phi = 0 there (one such end = torsion cantilever), solved by 1-dof
  // linear elements (nodal phi exact for this ODE); the torque diagram T(x)
  // then follows by statics, reusing sfdBmd (its V output).
  const twistMM=twistEnds(S);
  const twistSupportsMM=twistMM.map(t=>({pos:t.pos,type:'pinned',end:t.end}));
  let tors=null;
  const swE=selfWeightEccentricity(sec);
  const anyUserEcc = S.eccOn && S.loads.some(ld=>!ld.isSelfWeight && ld.type!=='moment' && Math.abs(ld.e||0)>1e-9);
  const anySelfWeightEcc = S.eccOn && Math.abs(swE)>1e-9 && selfWeightValue(sec)>0 &&
    [...ulsCombos,...slsCombos].some(cb=>Math.abs(cb.factors.G??0)>1e-12);
  const anyEcc = anyUserEcc || anySelfWeightEcc;
  const torsErr = (anyEcc && !(sec.J>0))? "the section torsional constant I_T is zero or undefined in the section data" : null;
  // Torque list of ONE combination (pattern-aware through comboLoadPieces):
  // point torques P*e [N.mm] and distributed torques w*e [N.mm/mm]; applied
  // moments carry no torque; the per-load z_g is a load height, not a torque.
  // Shared by the St Venant FE (torque diagram) and the warping-torsion FE.
  const mkT=(combo)=>{ const out=[]; comboLoadPieces(combo).forEach(p=>{
      const f=p.factor;
      const le=p.e; // this load's own shear-centre offset, mm
      if(p.type==='point') out.push({type:'point',pos:p.pos,P:p.P*f*1000*le}); // kN -> N, x e mm -> N.mm
      else if(p.type==='udl'||p.type==='trap') out.push({type:'udl',x1:p.x1,x2:p.x2,w1:p.w1*f*le,w2:p.w2*f*le});
    });
      const gF=combo.factors.G??0, sw=selfWeightValue(sec);
      if(Math.abs(gF)>1e-12 && Math.abs(swE)>1e-9 && sw>0){
        out.push({type:'udl',x1:0,x2:S.L*1000,w1:sw*gF*swE,w2:sw*gF*swE,isAutoSelfWeight:true});
      }
      return out; };
  if(anyEcc && !torsErr){
    const GIt=81000*sec.J*1e4; // N.mm2 (G = 81000 N/mm2 per SN003a / P385)
    const solveT=(tq)=>{
      const nodes=buildNodes(L,twistSupportsMM,tq,120);
      const n=nodes.length;
      const K=Array.from({length:n},()=>new Array(n).fill(0));
      const F=new Array(n).fill(0);
      for(let el=0;el<n-1;el++){ const Le=nodes[el+1]-nodes[el], k=GIt/Le;
        K[el][el]+=k; K[el][el+1]-=k; K[el+1][el]-=k; K[el+1][el+1]+=k; }
      const idx=new Map(nodes.map((x,i)=>[+x.toFixed(6),i]));
      tq.forEach(ld=>{ if(ld.type==='point'){ const i=idx.get(+(+ld.pos).toFixed(6)); if(i!=null) F[i]+=ld.P; }
        else if(ld.type==='udl'){ for(let el=0;el<n-1;el++){ const xa=nodes[el],xb=nodes[el+1];
          if(xb<=ld.x1+1e-9||xa>=ld.x2-1e-9) continue; const Le=xb-xa;
          const tv=x=>{ if(ld.x2===ld.x1) return ld.w1; const s=(x-ld.x1)/(ld.x2-ld.x1); return ld.w1+(ld.w2-ld.w1)*s; };
          const ta=tv(xa),tb=tv(xb);
          F[el]+=Le*(2*ta+tb)/6; F[el+1]+=Le*(ta+2*tb)/6; } } });
      const fixed=new Set(); twistSupportsMM.forEach(s=>{ const i=idx.get(+(+s.pos).toFixed(6)); if(i!=null) fixed.add(i); });
      if(!fixed.size) throw 'Torsion: no end restrains twist R_x (torsional mechanism).';
      const free=[]; for(let d2=0;d2<n;d2++) if(!fixed.has(d2)) free.push(d2);
      const phi=new Array(n).fill(0);
      if(free.length){ const Kff=free.map(r=>free.map(cc=>K[r][cc])), Ff=free.map(r=>F[r]);
        const df=linsolve(Kff,Ff); free.forEach((dof,j)=>phi[dof]=df[j]); }
      const R=new Array(n).fill(0);
      for(let i=0;i<n;i++){ let s2=0; for(let j=0;j<n;j++) s2+=K[i][j]*phi[j]; R[i]=s2-F[i]; }
      const reactions=twistSupportsMM.map(s=>({pos:s.pos,type:'pinned',V:R[idx.get(+(+s.pos).toFixed(6))]}));
      const fb=sfdBmd(L,twistSupportsMM,tq,reactions);
      let Tm=0,Tp=0; fb.V.forEach((v,i)=>{ if(Math.abs(v)>Math.abs(Tm)){Tm=v;Tp=fb.xs[i];} });
      let pm=0,pp=0; phi.forEach((v,i)=>{ if(Math.abs(v)>Math.abs(pm)){pm=v;pp=nodes[i];} });
      return {nodes,phi,xs:fb.xs,T:fb.V,Tmax:Tm,Tpos:Tp,phiMax:pm,phiPos:pp};
    };
    const uls=ulsCombos.map(cb=>({combo:cb,r:solveT(mkT(cb))}));
    let gT=uls[0]; uls.forEach(u=>{ if(Math.abs(u.r.Tmax)>Math.abs(gT.r.Tmax)) gT=u; });
    const sls=slsCombos.map(cb=>({combo:cb,r:solveT(mkT(cb))}));
    let gW=sls[0]; sls.forEach(u=>{ if(Math.abs(u.r.phiMax)>Math.abs(gW.r.phiMax)) gW=u; });
    tors={on:true, GIt,
      Tmax:Math.abs(gT.r.Tmax)/1e6, Tpos:gT.r.Tpos/1000, governT:gT.combo.label,
      diag:{xs:gT.r.xs.map(x=>x/1000), T:gT.r.T.map(t=>t/1e6)},
      uls,
      TmaxSLS:Math.abs(gW.r.Tmax)/1e6, phiMax:Math.abs(gW.r.phiMax), phiPos:gW.r.phiPos/1000, governTw:gW.combo.label};
  }

  // ---- P385 open-section torsion: Method B closed forms where they apply
  //      (single fork-fork span, full-span distributed and/or point torques,
  //      warping free at both ends), otherwise the general warping-torsion FE
  //      of js/checks/torsion-fe.js (19 Sep 2026 gap closure, group G4) ----
  let torsO=null;
  if(anyEcc && torsErr && !sec.isBox){ torsO={ok:false,reason:torsErr}; }
  if(anyEcc && !torsErr && !sec.isBox){
    const tp=sec.tp;
    const IT=((tp&&tp.IT)? tp.IT : sec.J)*1e4;                     // mm4, P385 App A preferred
    const IwO=(((tp&&tp.Iw!=null)? tp.Iw : sec.Iw)||0)*1e12;       // mm6
    const GItO=81000*IT;
    const aa=IwO>0? Math.sqrt(E*IwO/GItO) : 0; // use the same E as the bending analysis
    // P385 Appendix C: fork ends at both ends of the span - lateral translation
    // U_y and twist R_x held at BOTH ends, warping free at both
    const forkBoth = ends.every(e=>e.uy&&e.rx);
    const warpAny = ends.some(e=>e.rx&&e.warp);
    const mk385=(combo)=>{
      const list=[];
      for(const p of comboLoadPieces(combo)){
        if(p.type==='moment') continue;
        const f=p.factor, le=p.e;
        if(!f||Math.abs(le)<1e-9) continue;
        if(p.type==='point'){ list.push({kind:'point',alpha:p.pos/L,T:p.P*f*1000*le}); }
        else {
          const x1=p.x1, x2=p.x2;
          if(x1>1e-6 || Math.abs(x2-L)>1e-6) return {ok:false,reason:'partial-span eccentric distributed load (the P385 closed forms Cases 3/4/10 cover full-span distributed torque only)'};
          const w1=p.w1, w2=p.w2;
          const wu=Math.min(w1,w2), dv=w2-w1;
          if(Math.abs(wu)>1e-12) list.push({kind:'ud',T:wu*f*le*L});
          if(Math.abs(dv)>1e-12) list.push({kind:'lin',T:Math.abs(dv)/2*f*le*L*Math.sign(dv*1), mirror:dv<0});
        }
      }
      const gF=combo.factors.G??0, sw=selfWeightValue(sec);
      if(Math.abs(gF)>1e-12 && Math.abs(swE)>1e-9 && sw>0){
        list.push({kind:'ud',T:sw*gF*swE*L});
      }
      return {ok:true,list};
    };
    // Why the closed forms cannot be used (empty = they can): layout, a
    // warping-fixed support, or a partial-span distributed torque in any
    // analysed combination.
    const feReasons=[];
    if(!forkBoth) feReasons.push(twistMM.length===1? 'twist restrained at End '+twistMM[0].end+' only (torsion cantilever)' : 'the ends are not both fork ends (U_y + R_x at both ends)');
    if(warpAny) feReasons.push('warping-fixed end');
    if(!feReasons.length){
      for(const cb of [...ulsCombos,...slsCombos]){ const m=mk385(cb); if(!m.ok){ feReasons.push(m.reason); break; } }
    }
    if(!(IT>0)) torsO={ok:false,reason:'the torsional constant I_T is zero or undefined'};
    else if(!(IwO>0)||!(aa>0)||!isFinite(aa)) torsO={ok:false,reason:'no warping constant available for this section'};
    else if(!feReasons.length){
      // SCI P385 Appendix C closed forms, superposed per combination
      const sols=ulsResults.map(res=>({combo:res.combo, fb:res.fb, sol:p385Solve(L,aa,GItO,mk385(res.combo).list)}));
      let slsSol=null;
      for(const cb of slsCombos){
        const s2=p385Solve(L,aa,GItO,mk385(cb).list);
        let pm=0,pp=0; s2.phi.forEach((v,i)=>{ if(Math.abs(v)>Math.abs(pm)){pm=v;pp=s2.xs[i];} });
        if(!slsSol||Math.abs(pm)>Math.abs(slsSol.phiMax)) slsSol={combo:cb,phiMax:pm,phiPos:pp};
      }
      torsO={ok:true,method:'closed',methodLabel:'SCI P385 App C closed forms (Cases 3/4/10)',
        bcText:'fork ends at x = 0 and x = L (&phi; = 0, warping free)',feReasons:[],
        aa,X:L/aa,IT,Iw:IwO,GIt:GItO,sols,sls:slsSol};
    } else {
      // General warping-torsion FE: E I_w phi'''' - G I_T phi'' = m_t(x), phi = 0
      // at every end with R_x held, phi' = 0 where warping is restrained, a free
      // end natural; mesh doubled once for the error.
      const feSup=twistMM.map(t=>({pos:t.pos,warpFix:t.warpFix,end:t.end}));
      const EIwO=E*IwO;
      const feStats={solved:0,cached:0};
      const feSolve=(cb)=>warpingTorsionFEMemo({L,EIw:EIwO,GIt:GItO,supports:feSup,torques:mkT(cb)},feStats);
      const sols=ulsResults.map(res=>({combo:res.combo, fb:res.fb, sol:feSolve(res.combo)}));
      let slsSol=null, meshError=0, nElem=0, nElemCoarse=0;
      sols.forEach(se=>{ meshError=Math.max(meshError,se.sol.meshError); nElem=Math.max(nElem,se.sol.nElem); nElemCoarse=Math.max(nElemCoarse,se.sol.nElemCoarse); });
      for(const cb of slsCombos){
        const s2=feSolve(cb);
        meshError=Math.max(meshError,s2.meshError); nElem=Math.max(nElem,s2.nElem); nElemCoarse=Math.max(nElemCoarse,s2.nElemCoarse);
        let pm=0,pp=0; s2.phi.forEach((v,i)=>{ if(Math.abs(v)>Math.abs(pm)){pm=v;pp=s2.xs[i];} });
        if(!slsSol||Math.abs(pm)>Math.abs(slsSol.phiMax)) slsSol={combo:cb,phiMax:pm,phiPos:pp};
      }
      torsO={ok:true,method:'fe',methodLabel:'warping-torsion FE ('+nElem+' elements)',nElem,nElemCoarse,meshError,
        converged:meshError<=TORSION_FE_MESH_BLOCK,meshBlock:TORSION_FE_MESH_BLOCK,nSolves:feStats.solved,nCached:feStats.cached,
        bcText:torsionFeBcText(feSup,L),feReasons,
        aa,X:L/aa,IT,Iw:IwO,GIt:GItO,sols,sls:slsSol};
    }
  }
  return {sec,py,E,L,Ix,tors,torsO,torsErr,swPerM:sec.mass*9.81/1000,ulsResults,
    Vmax:Vmax/1000, Mmax:Mmax/1e6, Mpos:Mpos/1000,
    Mq:Mq/1e6, Mh:Mh/1e6, Mq3:Mq3/1e6, M24:M24/1e6, M0end:M0end/1e6, MLend:MLend/1e6,
    dmax, dpos:dpos/1000, deflection:governD.deflection, deflSegments, divisorCant, deflAbs,
    diag:{xs:xs.map(x=>x/1000), V:Venv.map(v=>v/1000), M:Menv.map(m=>m/1e6),
          dx:governD.r.nodes.map(x=>x/1000), dw:governD.r.w},
    reactions, ulsResults, ulsCompanions, slsResults, governV, governM, governD,
    ends, supports:endsToSupports(S), cant, stability, companionNote, ulsCombos, slsCombos, companionCombos,
    uplift:{list:uplift, supports:upliftSupports, any:upliftSupports.length>0, nCombos:ulsResults.length+ulsCompanions.length+slsResults.length}};
}

/* ---- Hold-down check (19 Sep 2026 gap closure, item 1.2), pure ----
   One row per end that lifts in any combination (n = end number). An end
   lifting in a ULS combination (the gamma_G,inf companions included) is
   blocking (unsupported) unless its "hold-down provided" box is ticked, in
   which case it is an advisory carrying the design force. An end that lifts
   ONLY in SLS combinations (a deflection case is not a ULS equilibrium
   state; since 20 Sep 2026 the default SLS combination is 1.0G + 1.0Q, so an
   SLS-only lift needs a user combination without G, e.g. the NA 2.23
   "variable actions only" case) is reported as an advisory naming the
   combination and the force.
   Returns {rows, unsupported, advisory}; rows = [{n, pos (mm), R (kN,
   negative), combo, sls, RUls, comboUls, RSls, comboSls, holdDown, level:
   'uls'|'sls', blocking, msg}]. */
function holdDownCheck(a){
  const out={rows:[],unsupported:[],advisory:[]};
  const up=a.uplift&&a.uplift.supports||[];
  up.forEach(u=>{
    const kN=v=>(Math.abs(v)).toFixed(2);
    const where='at End '+u.n+' (x = '+(u.pos/1000).toFixed(2).replace(/\.?0+$/,'')+' m)';
    let msg, blocking=false, level;
    if(u.RUls!=null){
      level='uls';
      const force=(withWhere)=>'R = &minus;'+kN(u.RUls)+' kN'+(withWhere? ' '+where : '')+' (combination '+u.comboUls+')'+(u.RSls!=null? '; SLS uplift &minus;'+kN(u.RSls)+' kN ('+u.comboSls+')' : '');
      if(u.holdDown) msg='Hold-down provided '+where+': design the hold-down for '+force(false)+'. Reaction taken as tension at the support; the connection and the supporting structure are not designed here.';
      else { blocking=true; msg='Hold-down required: '+force(true)+'. The end cannot resist uplift as modelled; tick "hold-down provided" for this end once a holding-down connection is designed for this force, or revise the restraints / loading (EN 1990 2.4.4 EQU; the &gamma;<sub>G,inf</sub> companions - G at 1.0 (STR set B) and 0.9 (EQU set A) with the entered variable factors - are included).'; }
    } else {
      level='sls';
      msg='Hold-down check (SLS only) '+where+': the SLS combination '+u.comboSls+' lifts this end by R = &minus;'+kN(u.RSls)+' kN; no ULS combination lifts it, including the &gamma;<sub>G,inf</sub> companions with G at 1.0 (STR set B) and 0.9 (EQU set A, EN 1990 Table A1.2(A)) and the entered variable factors'+(u.holdDown? '; hold-down provided' : '')+'. A deflection (SLS) combination is not a ULS equilibrium state, so this does not block PASS.';   // wording only (20 Sep 2026: the SLS default now carries G)
    }
    const row=Object.assign({},u,{msg,level,blocking});
    out.rows.push(row);
    if(blocking) out.unsupported.push(msg); else out.advisory.push(msg);
  });
  return out;
}

/* Standard-specific check engines live in js/checks/. */
function analysisForCombination(a,res){
  const fb=res.fb, L=a.L;
  let M24=0; fb.xs.forEach((x,i)=>{ if(x>=L/4&&x<=3*L/4) M24=Math.max(M24,Math.abs(fb.M[i])); });
  return Object.assign({},a,{ulsResults:[res],governM:res,governV:res,
    Mmax:res.Mmax/1e6,Vmax:res.Vmax/1000,Mpos:res.Mpos/1000,
    Mq:interpAt(fb.xs,fb.M,L/4)/1e6,Mh:interpAt(fb.xs,fb.M,L/2)/1e6,
    Mq3:interpAt(fb.xs,fb.M,3*L/4)/1e6,M24:M24/1e6,
    M0end:interpAt(fb.xs,fb.M,1e-4)/1e6,MLend:interpAt(fb.xs,fb.M,L-1e-4)/1e6});
}
function checks(a){
  if(S.code==='EC3') return (S.restraint||'full')==='full'? checksEC3Restrained(a) : checksEC3UnrestrainedSCI(a);
  // mLT and mx depend on each diagram; a lower peak moment can govern.
  const results=a.ulsResults.map(res=>({res,c:checksBS5950(analysisForCombination(a,res))}));
  const c=checksBS5950(a); // retain envelope quantities for the detailed report
  c.utils=c.utils.map((u,i)=>{
    const combo=u.name==='Deflection'?a.governD.combo.label:u.name.startsWith('Shear')?a.governV.combo.label:a.ulsResults.length===1?a.governM.combo.label:'envelope';
    let worst={...u,combo};
    results.forEach(r=>{ const v=r.c.utils[i]; if(v&&v.val>worst.val) worst={...v,combo:r.res.combo.label}; });
    return worst;
  });
  c.unsupported=[...new Set(c.unsupported.concat(...results.map(r=>r.c.unsupported)))];
  // shared analysis-level checks (uplift / hold-down, stability notes, the gamma_G,inf companion note)
  const hd=holdDownCheck(a);
  c.holdDown=hd;
  c.unsupported=c.unsupported.concat(hd.unsupported);
  c.advisory=(c.advisory||[]).concat(hd.advisory);
  if(a.stability&&a.stability.notes) c.advisory.push(...a.stability.notes);
  if(a.companionNote) c.advisory.push(a.companionNote);
  c.gov=c.utils.reduce((p,u)=>u.val>p.val?u:p);
  c.pass=c.unsupported.length===0&&c.utils.every(u=>Number.isFinite(u.val)&&u.val>=0&&u.val<=1.0001);
  c.combinationChecks=results.map(r=>({combo:r.res.combo.label,utils:r.c.utils}));
  return c;
}



/* ==== beam-v03 module 10 ==== */
/* ---------------------------------------------------------------------------
   BS 5950-1:2000 LTB effective length from the end flags (pure; reads
   st.ends, st.leFactor, st.destab; 19 Sep 2026 review finding: a released
   U_y / R_x passed on L_E = 1.0 L). Returns {LE (mm), K, basis (HTML,
   [verify]), row, msg (HTML, blocking) | null}.
   - An entered L_E/L factor: K = factor x (1.2 when destabilising), as before
     (the user's own Table 13 / 14 reading).
   - Fork ends at both ends (U_y + R_x held): Table 13, compression flange
     laterally restrained, nominal torsional restraint - "both flanges free to
     rotate on plan" (R_z free at both ends) 1.0 L normal / 1.2 L
     destabilising; "both flanges fully restrained against rotation on plan"
     (R_z held at both ends) 0.7 L / 0.85 L; R_z at one end only: Table 13 has
     no one-end row, the 1.0 L / 1.2 L row is taken (conservative).
   - Cantilever (in-plane fixed root - free tip) whose root holds U_y, R_z and
     R_x: Table 14 (d) "restrained laterally, torsionally and against
     rotation on plan" (the root warping flag is not distinguished by the
     table), tip: free 0.8 L / 1.4 L; torsional restraint (R_x) 0.6 L / 0.6 L;
     lateral and torsional restraint (U_y + R_x) 0.5 L / 0.5 L; a tip holding
     U_y only is taken as free (the Table 14 "lateral restraint to top flange"
     row is not claimed for a shear-centre restraint: conservative).
   - Any other end flag set (an end releasing U_y or R_x, a lateral
     cantilever, a cantilever root releasing R_z): no tabulated row - blocked
     until an L_E/L factor is entered.
   The destabilising column of the table is used directly (not the x1.2 device
   of the entered-factor path).
   --------------------------------------------------------------------------- */
function bs5950LtbLength(st,L){
  st=st||S;
  const [e1,e2]=endsList(st);
  const destab=!!st.destab;
  const sym={uy:'U<sub>y</sub>',rx:'R<sub>x</sub>'};
  if(st.leFactor!=null&&st.leFactor!==''&&Number.isFinite(+st.leFactor)){
    const K=(+st.leFactor)*(destab?1.2:1);
    return {LE:K*L, K, row:'user', msg:null,
      basis:'entered L<sub>E</sub>/L = '+(+st.leFactor).toFixed(2)+(destab? ' &times; 1.2 (destabilising)' : '')+' (the end fixities would give: '+bs5950LtbLength(Object.assign({},st,{leFactor:null}),L).basis+')'};
  }
  const cantIP=!!(e1.uz&&e1.ry&&!e2.uz&&!e2.ry);   // in-plane cantilever: Table 14, never Table 13
  if(cantIP && e1.uy&&e1.rz&&e1.rx){
    const tip= (e2.uy&&e2.rx)? 'lateral and torsional restraint (U<sub>y</sub> + R<sub>x</sub>)' : e2.rx? 'torsional restraint (R<sub>x</sub>)' : 'free'+(e2.uy? ' (U<sub>y</sub> at the tip not claimed as the top-flange row: conservative)' : '');
    const K= (e2.uy&&e2.rx)? 0.5 : e2.rx? 0.6 : (destab? 1.4 : 0.8);
    const row='Table 14 (d): root restrained laterally, torsionally and against rotation on plan; tip '+tip;
    return {LE:K*L, K, row, msg:null, basis:K.toFixed(2)+' L from the end flags - cantilever, '+row+', '+(destab? 'destabilising' : 'normal')+' loading (the root warping flag is not distinguished by Table 14) [verify]'};
  }
  const forkBoth=!!(e1.uy&&e1.rx&&e2.uy&&e2.rx);
  if(!cantIP && forkBoth){
    const rzBoth=!!(e1.rz&&e2.rz), rzOne=!rzBoth&&!!(e1.rz||e2.rz);
    const K= rzBoth? (destab? 0.85 : 0.7) : (destab? 1.2 : 1.0);
    const row= rzBoth? 'Table 13: both flanges fully restrained against rotation on plan (R<sub>z</sub> held at both ends)' : 'Table 13: both flanges free to rotate on plan'+(rzOne? ' (R<sub>z</sub> held at one end only: no one-end row, taken as free to rotate - conservative)' : '');
    return {LE:K*L, K, row, msg:null, basis:K.toFixed(2)+' L from the end flags - fork ends (U<sub>y</sub> + R<sub>x</sub>) at both ends, compression flange laterally restrained with nominal torsional restraint, '+row+', '+(destab? 'destabilising' : 'normal')+' loading [verify]'};
  }
  const rel=[e1,e2].map(e=>{ const r=['uy','rx'].filter(k=>!e[k]).map(k=>sym[k]); return r.length? 'End '+e.n+' releases '+r.join(' and ') : null; }).filter(Boolean);
  const why= cantIP? 'the cantilever root (End 1) does not hold U<sub>y</sub>, R<sub>z</sub> and R<sub>x</sub> (Table 14 (d))' : rel.join('; ');
  return {LE:1.0*L, K:1.0, row:'none', basis:'no tabulated row for these end flags (1.0 L printed, not accepted)',
    msg:'BS 5950 lateral-torsional buckling: the end flags have no Table 13 / Table 14 row - '+why+'. L<sub>E</sub> cannot be taken as 1.0 L for a released end; enter the L<sub>E</sub>/L factor for these end conditions (Table 13 for a beam with both ends held, Table 14 for a cantilever) under Axial &amp; lateral-torsional buckling. PASS is blocked.'};
}
function checksBS5950(a){
  const sec=a.sec, py=a.py, E=a.E;
  const unsupported=[];
  const advisory=[];
  const eps=Math.sqrt(275/py);
  const Ag=sec.A*1e2;
  const F=S.axial||0, Fc0=Math.max(F,0);
  if(Math.abs(S.Mz||0)>1e-9) unsupported.push('BS 5950 minor-axis bending and biaxial interaction are not implemented. The entered minor-axis moment cannot be ignored; PASS is blocked.');
  // BS 5950-1:2000 Figure 5 flange outstand: rolled I/H  b = B/2; channel  b = B
  // (the FULL flange width). The section tables store the EC3 ratio
  // c/t = ((B - tw - 2r)/2)/tf, which is SMALLER and would misclassify near the
  // Table 11 limits, so the BS ratio is rebuilt from the raw geometry here.
  // The web d/t is the same dimension in both codes (d = depth between fillets),
  // and the box tables already store the Table 12 flats (B-3t HF / B-5t CF per
  // note a), so those pass through unchanged.
  // Classification sees the coexistent axial compression: Table 11 "generally"
  // web rows with r1/r2 (cl 3.5.5) tighten the limits below 80/100/120e.
  const bT_BS = sec.isBox? sec.bT : ((sec.kind==='channel'? sec.B : sec.B/2)/sec.tf);
  const cl = sec.isBox? classifyBox(sec.bT,sec.dt,eps,sec.boxType)
    : classify(bT_BS,sec.dt,eps,sec.kind,{Fc:Fc0*1000, dwt:(sec.d||0)*sec.tw, Ag:Ag, py:py});
  const clsName=["","Plastic","Compact","Semi-compact","Slender"][cl.cls];
  if(cl.cls>=4) unsupported.push("BS 5950 slender section"+(cl.webCase==='bending+compression'?" (web classified for combined bending + compression, Table 11 r1/r2 rows)":"")+": effective-section design is required; gross Zx is not accepted for PASS.");
  if(sec.isBox && Fc0>0) advisory.push("Box-section web classification uses the pure-bending Table 12 limits; the r1-based rows for coexistent axial compression are not applied (axial + bending on boxes is blocked from PASS separately when significant).");
  advisory.push("Web bearing and buckling at supports and under concentrated loads (cl 4.5.2 / 4.5.3) are outside this calculator's scope - verify separately wherever a reaction or point load is applied to an unstiffened web.");
  // shear   BS 5950 cl 4.2.3: open section Av=tw D; box section Av=A D/(D+B)
  // (for a square SHS, D=B, so this reduces exactly to A/2 as before).
  const Av = sec.isBox? sec.A*100*sec.D/(sec.D+sec.B) : sec.tw*sec.D; // sec.A is cm  for box ? mm ; PFC/UB/UC path unaffected
  const Pv=0.6*py*Av/1000, Fv=Math.abs(a.Vmax);
  const lowShear=Fv<=0.6*Pv;
  const shearBuckle=sec.dt>70*eps;
  if(shearBuckle) unsupported.push("BS 5950 shear buckling check is required and is not implemented in this calculator.");
  // moment capacity
  const Zx=sec.Zx*1e3, Sx=sec.Sx*1e3;
  let Mcx;
  if(cl.cls<=2) Mcx=Math.min(py*Sx,1.2*py*Zx)/1e6;
  else Mcx=py*Zx/1e6;
  let hsNote=null;
  if(!lowShear){
    if(sec.isBox){
      unsupported.push("BS 5950 high-shear moment reduction for box/RHS/SHS sections is section-specific and is not implemented exactly.");
    } else if(cl.cls<=2){
      const rho=Math.pow(2*Fv/Pv-1,2);
      const Sv=Av*Av/(4*sec.tw); // mm3, web shear area plastic modulus for open sections
      Mcx=Math.min(py*(Sx-rho*Sv),1.2*py*Zx)/1e6;
      hsNote=`high-shear reduction applied with S<sub>v</sub>=A<sub>v</sub> /(4t<sub>w</sub>)`;
    } else {
      // cl 4.2.5.3, class 3 semi-compact: Mc = py(Z - rho*Sv/1.5)
      const rho=Math.pow(2*Fv/Pv-1,2);
      const Sv=Av*Av/(4*sec.tw);
      Mcx=py*Math.max(Zx-rho*Sv/1.5,0)/1e6;
      hsNote=`high-shear reduction (semi-compact, cl 4.2.5.3): M<sub>c</sub> = p<sub>y</sub>(Z<sub>x</sub> &minus; &rho;S<sub>v</sub>/1.5)`;
    }
  }
  // effective area & axial
  const Anet=(S.anet!=null? S.anet*1e2 : Ag);
  const Ke=(S.Ke!=null? S.Ke : (KeByGrade[S.grade]||1.2));
  const Ae=Math.min(Ke*Anet,Ag);
  const Pz=Ae*py/1000;
  const n=Pz>0? Math.abs(F)/Pz : 0;
  // Reduced modulus for coexistent axial (cl 4.2.5 / 4.8.3.2). Class 1/2 uses
  // the plastic web-block formula S_r = S - (n Ag)^2/(4 t_w), valid only while
  // the shifted plastic neutral axis stays IN THE WEB (n Ag <= t_w d) - beyond
  // that the formula understates the loss, so PASS is blocked. Class 3 uses the
  // ELASTIC basis Z_x(1 - n) (linear interaction, cl 4.8.3.2) - previously the
  // plastic S_x was used for every class, overstating semi-compact capacity by
  // the shape factor (~12-15% on a UB). Only the covered I/H path feeds PASS;
  // box/channel axial+bending stays blocked.
  let Srx = cl.cls<=2? Sx : Zx;
  const tEff = sec.isBox? 2*sec.tw : sec.tw;
  if(n>0.02){
    if(sec.isBox || sec.kind==='channel') unsupported.push("BS 5950 axial-load plus major-axis bending requires section-family-specific reduced modulus data; this calculator does not use an adapted approximate formula for PASS.");
    else if(cl.cls<=2){
      const dn=n*Ag/sec.tw; // depth of web carrying the axial block, mm
      if(dn>(sec.d||0)+1e-6) unsupported.push("Axial ratio n = "+n.toFixed(2)+" pushes the plastic neutral axis out of the web (n&middot;A<sub>g</sub>/t<sub>w</sub> &gt; d): the web-block reduced-modulus formula is no longer valid and the flange-region expression is not implemented; PASS is blocked.");
      Srx=Math.max(0, Sx-(Ag*Ag*n*n)/(4*tEff));
    } else {
      Srx=Math.max(0, Zx*(1-n));
    }
  }
  const Mrx=Math.max(0,Math.min(Mcx,(cl.cls<=2? Math.min(py*Srx, 1.2*py*Zx) : py*Srx)/1e6));
  if(Math.abs(F)>1e-9&&!lowShear) unsupported.push('BS 5950 combined high shear, axial force and bending requires the web interaction check in clause 4.8; this is not implemented. PASS is blocked.');
  const Mx=Math.abs(a.Mmax);
  // zero/negative reduced capacity with a coexistent moment must read as a
  // failure, not util = 0 (n >= 1 zeroes the Class 3 elastic basis)
  const localUtil=Mrx>1e-9? Mx/Mrx : (Mx>1e-9? 99 : 0);
  // m-factors. Table 18 note: mLT = 1.0 for cantilevers AND for members with
  // DESTABILISING loading conditions (previously the destabilising switch only
  // lengthened LE and kept mLT < 1, which Table 18 does not permit).
  const isCant=isCantilever(S);
  const mf=mFactors(a.Mq,a.Mh,a.Mq3,a.Mmax,a.M24);
  let mLT=(isCant||S.destab)?1:mf.mLT, mx=mf.mx;
  if(S.mLTo!=null) mLT=S.mLTo;
  if(S.mxo!=null) mx=S.mxo;
  // LTB   BS 5950 box-section path. SHS naturally returns very low ?LT because
  // Ix Iy; RHS uses the closed-section ?LT expression rather than a rough Table
  // 15 screen.
  // LTB effective length: the entered L_E/L factor (x1.2 destabilising) or, when
  // blank, Table 13 / 14 from the end flags (bs5950LtbLength; blocked where no
  // row exists - 19 Sep 2026 review)
  const leBS=bs5950LtbLength(S,a.L);
  const LE=leBS.LE;
  if(leBS.msg) unsupported.push(leBS.msg);
  const ry=sec.ry*10;
  let lam=null,v=null,betaW=null,lamLT=null,pb=null,lamL0=null,Mb,ltbUtil,rhsFlag=false,phiB=null,gammaPrime=null;
  if(sec.isBox){
    lam=LE/ry; lamL0=0.4*Math.sqrt(Math.PI*Math.PI*E/py);
    betaW=cl.cls<=2?1:Zx/Sx;
    const Ixmm=sec.Ix*1e4, Iymm=sec.Iy*1e4, Jmm=sec.J*1e4;
    gammaPrime=Math.max(0,(1-Iymm/Ixmm)*(1-Jmm/(2.6*Ixmm)));
    phiB=Math.sqrt(Math.max((Sx*Sx*gammaPrime)/(Ag*Jmm),0));
    lamLT=2.25*Math.sqrt(Math.max(phiB*lam*betaW,0));
    ({pb,lamL0}=pbFunc(lamLT,py,E));
    Mb=Math.min((cl.cls<=2? pb*Sx : pb*Zx)/1e6, Mcx);
    mLT=1; ltbUtil=Mb>0? Mx/Mb : 0;
  } else {
    lam=LE/ry;
    v=1/Math.pow(1+0.05*Math.pow(lam/sec.x,2),0.25);
    betaW=cl.cls<=2?1:Zx/Sx;
    lamLT=sec.u*v*lam*Math.sqrt(betaW);
    ({pb,lamL0}=pbFunc(lamLT,py,E));
    Mb=(cl.cls<=2? pb*Sx : pb*Zx)/1e6;
    // cl 4.3.6.2: the buckling resistance check is mLT*Mx <= Mb (the companion
    // Mx <= Mcx is carried by the bending utilisation above). Checking plain
    // Mx <= Mb was over-conservative by 1/mLT; mLT = 1 for cantilevers and
    // destabilising loads, so those cases are unchanged.
    ltbUtil=Mb>0? mLT*Mx/Mb : 0;
  }
  // strut (axial term)   rx=ry for a square SHS, so Pc=Pcy automatically there;
  // for UB the major/minor axis curves genuinely differ (Table 23).
  const autoRob = defaultRobertson(S.family,sec.boxType,sec.tf);
  const a_robX = S.robX!=null? S.robX : autoRob.x;
  const a_robY = S.robY!=null? S.robY : autoRob.y;
  const rx=sec.rx*10;
  // strut lengths: the entered L_E/L factor, else from the end fixities
  // (lcrDefaults, Table 22 style, x-x from U_z / R_y and y-y from U_y / R_z);
  // the destabilising x1.2 is kept on the strut length as before on this path
  const lcr=lcrDefaults(S);
  const Fc=Math.max(F,0);
  if(Fc>0 && (lcr.Ky==null || lcr.Kz==null)) throw 'Strut buckling: the end fixities form a mechanism ('+(lcr.Ky==null? lcr.basisY : lcr.basisZ)+'); F_c cannot be carried.';
  const Kx=lcr.Ky!=null? lcr.Ky : 1.0, Ky=lcr.Kz!=null? lcr.Kz : 1.0;
  const LcrX=Kx*(S.destab?1.2:1)*a.L, LcrY=Ky*(S.destab?1.2:1)*a.L, lcrBasis=lcr.basis;
  const pcx=pcFunc(LcrX/rx,py,a_robX,E), pcy=pcFunc(LcrY/ry,py,a_robY,E);
  const Pc=Ag*pcx/1000, Pcy=Ag*pcy/1000;
  if(Fc>0 && sec.kind==='channel') unsupported.push("BS 5950 PFC/channel compression must use the UK channel strut approach/Table 25 or verified section tables data; the previous generic Robertson placeholder is not accepted.");
  const pyZx=py*Zx/1e6;
  const u1=Fc/Pc + mx*Mx/pyZx;
  const u2=Fc/Pcy + mLT*Mx/Mb;
  // deflection - governing segment with its own limit (span/divisor, L/divisorCant
  // for a cantilever segment, capped by the optional absolute limit), from analyse()
  const span=a.deflection?a.deflection.span:a.L;
  const divisor=(a.deflection&&a.deflection.divisor!=null)? a.deflection.divisor : S.divisor;
  const dlimit=(a.deflection&&a.deflection.limit!=null)? a.deflection.limit : span/divisor;
  const dmax=Math.abs(a.deflection?a.deflection.dmax:a.dmax), defOk=dmax<=dlimit;
  const deflCant=!!(a.deflection&&a.deflection.cant), deflAbsGoverns=!!(a.deflection&&a.deflection.absGoverns);

  if(S.eccOn && S.loads.some(ld=>Math.abs(ld.e||0)>1e-9)) unsupported.push("Load eccentricity / torsion design is implemented for the EC3 code path only; switch Design code to EC3.");
  const utils=[
    {name:"Shear  Fv/Pv",val:Fv/Pv},
    {name:"Bending  Mx/Mrx",val:localUtil},
    {name: sec.isBox? "Mx/Mcx":"LTB  mLT.Mx/Mb", val:ltbUtil},
    {name:"Buckling (in-plane)",val:u1},
    {name:"Buckling (LTB interaction)",val:u2},
    {name:"Deflection",val:dmax/dlimit},
  ];
  if(Math.abs(F)>1e-9) utils.push({name:F<0?'Tension  Ft/Pt':'Compression cross-section  Fc/(Ag.py)',val:Math.abs(F)/(F<0?Pz:Ag*py/1000)});
  let gov=utils[0]; utils.forEach(u=>{ if(u.val>gov.val) gov=u; });
  const pass=unsupported.length===0 && utils.every(u=>u.val<=1.0001);

  return {eps,cl,clsName,unsupported,advisory,bTBS:bT_BS,Av,Pv,Fv,lowShear,shearBuckle,Mcx,hsNote,Zx,Sx,rhsFlag,
    Ag,Anet,Ke,Ae,Pz,F,n,Srx,Mrx,Mx,localUtil,isCant,mLT,mx,mf,
    LE,leK:leBS.K,leBasis:leBS.basis,leRow:leBS.row,leMsg:leBS.msg,LcrX,LcrY,Kx,Ky,lcrBasis,lam,v,betaW,phiB,gammaPrime,lamLT,pb,lamL0,Mb,ltbUtil,a_robX,a_robY,pcx,pcy,Pc,Pcy,Fc,pyZx,u1,u2,
    span,divisor,dlimit,dmax,defOk,deflCant,deflAbsGoverns,utils,gov,pass};
}



/* ==== beam-v03 module 11 ==== */
function checksEC3(a){
  const sec=a.sec, fy=a.py, E=a.E; // fy reuses the same Table 9 / EN 10025-2 thickness bands
  const unsupported=[];
  const eps=epsEC3(fy);
  const gM0=1.0, gM1=1.0; // UK NA
  const cl=classifyEC3(sec,eps);
  const clsName=["","Class 1","Class 2","Class 3","Class 4"][cl.cls];
  if(cl.cls>=4) unsupported.push("EC3 Class 4 section: effective-section properties per EN 1993-1-5 are required; gross Wel is not accepted for PASS.");

  const Av=avEC3(sec);
  const VplRd=Av*fy/(Math.sqrt(3)*gM0)/1000; // kN
  const Fv=Math.abs(a.Vmax);
  const lowShear=Fv<=0.5*VplRd; // EC3 cl 6.2.8: high-shear threshold is 0.5 Vpl,Rd (not 0.6 like BS5950)
  const shearBuckle = !sec.isBox && (sec.dt/eps) > 72; // cl 6.2.6(6): web shear-buckling check needed if hw/t>72e/eta (eta~1.0 conservative)
  if(shearBuckle) unsupported.push("EC3 web shear-buckling check per EN 1993-1-5 is required and is not implemented in this calculator.");

  const Zx=sec.Zx*1e3, Sx=sec.Sx*1e3;
  let McRd = cl.cls<=2? Sx*fy/gM0/1e6 : Zx*fy/gM0/1e6;
  let hsNote=null;
  if(!lowShear){
    if(cl.cls<=2 && sec.kind==='I'){
      const rho=Math.pow(2*Fv/VplRd-1,2);
      const Sv=Av*Av/(4*sec.tw);
      McRd=Math.min(McRd,(Sx-rho*Sv)*fy/gM0/1e6);
      hsNote='high-shear reduction applied with A<sub>v</sub> /(4t<sub>w</sub>) for a rolled I/H section';
    } else {
      unsupported.push("Exact EC3 high-shear reduced moment resistance for this section family/class is not implemented.");
    }
  }

  const Ag=sec.A*1e2; // mm2 (gross; bolt-hole deduction not modelled here)
  const NRd=Ag*fy/gM0/1000; // kN
  const F=S.axial||0;
  const n=NRd>0? Math.abs(F)/NRd : 0;
  let MNRd=McRd, awNote=false;
  if(n>0.01 && cl.cls<=2){
    if(sec.kind==='channel'){
      unsupported.push("EC3 reduced moment resistance under axial force is not defined here for PFC/channel sections; use verified section tables/software data.");
    } else {
      const aw=Math.min(Math.max((Ag-2*sec.B*sec.tf)/Ag,0),0.5);
      MNRd=Math.min(McRd*(1-n)/(1-0.5*aw), McRd);
      awNote=true;
    }
  }
  const Mx=Math.abs(a.Mmax);
  const localUtil=MNRd>0? Mx/MNRd : 0;

  const isCant=isCantilever(S);
  const c1r = S.C1o!=null? {C1:S.C1o, method:'user override'} : computeC1(a,isCant);
  const C1=c1r.C1;

  const LE=ltbLeFactor()*(S.destab?1.2:1)*a.L;
  let MbRd,ltbUtil,lamLT=null,chiLT=null,kc=null,fmod=null,curveInfo=null,Mcr=null,chiLTmod=null;
  if(sec.isBox){
    MbRd=McRd; ltbUtil=MbRd>0? Mx/MbRd:0;
  } else {
    if(sec.kind==='channel' && Mx>1e-9) unsupported.push("EC3 PFC/channel LTB requires a channel-specific Mcr calculation including exact warping/shear-centre/load-position terms; the previous approximate Mcr path is not accepted for PASS.");
    Mcr=mcrEC3(sec,LE,E,fy,C1)/1e6; // kN.m
    const Wy=(cl.cls<=2? Sx : Zx);
    lamLT=Math.sqrt(Wy*fy/(Mcr*1e6));
    curveInfo=ltbCurveEC3(sec);
    chiLT=chiLTEC3(lamLT,curveInfo.alphaLT,0.4,0.75);
    kc=1/Math.sqrt(Math.max(C1,1e-6));
    fmod=Math.min(1-0.5*(1-kc)*(1-2*Math.pow(lamLT-0.8,2)),1.0);
    chiLTmod=Math.min(chiLT/fmod,1.0);
    MbRd=Math.min(chiLTmod*Wy*fy/gM1/1e6, McRd);
    ltbUtil=MbRd>0? Mx/MbRd:0;
  }

  const ry=sec.ry*10, rx=sec.rx*10;
  const lam1=Math.PI*Math.sqrt(E/fy);
  const curveY=strutCurveEC3(sec,'y'), curveZ=strutCurveEC3(sec,'z');
  const lamBarY=(LE/rx)/lam1, lamBarZ=(LE/ry)/lam1;
  const chiY=chiStrutEC3(lamBarY,curveY.alpha), chiZ=chiStrutEC3(lamBarZ,curveZ.alpha);
  const NbRdY=chiY*Ag*fy/gM1/1000, NbRdZ=chiZ*Ag*fy/gM1/1000;
  const Fc=Math.max(F,0);

  // Member-buckling interaction (cl 6.3.3) via Annex B, Method 2   the method
  // recommended by steelconstruction.info as "the simpler approach for manual
  // calculations" (Annex A is the alternative; both are permitted by the UK NA).
  // Table B.1 applies to members NOT susceptible to torsional deformation
  // (CHS/SHS/RHS); Table B.2 applies to members that ARE susceptible (I, H,
  // channel). The two tables share the same kyy/kzz, but differ in kyz/kzy.
  const n_ratio = Fc/Math.max(Ag*fy/gM0/1000,1e-9); // NEd/NRd
  const psiInteraction = a.M0end!==0? Math.max(-1,Math.min(1,a.MLend/(a.M0end||1e-9))) : 0;
  const endMomentOnly = Math.max(Math.abs(a.M0end),Math.abs(a.MLend)) > 0.98*Math.max(Math.abs(a.Mmax),1e-9);
  const Cm = endMomentOnly ? Math.max(0.6+0.4*psiInteraction, 0.4) : 1.0;
  const cmMethod = endMomentOnly ? 'linear end-moment diagram' : 'general transverse-load moment diagram; conservative Cm=1.0';
  const Cmy=Cm, Cmz=Cm, CmLT=Cm;

  const ny=Fc/Math.max(NbRdY,1e-9), nz=Fc/Math.max(NbRdZ,1e-9);
  let kyy,kzz,kyz,kzy;
  if(cl.cls<=2){
    kyy=Math.min(Cmy*(1+(lamBarY-0.2)*ny), Cmy*(1+0.8*ny));
    kzz=Math.min(Cmz*(1+(2*lamBarZ-0.6)*nz), Cmz*(1+1.4*nz));
  } else {
    kyy=Math.min(Cmy*(1+0.6*lamBarY*ny), Cmy*(1+0.6*ny));
    kzz=Math.min(Cmz*(1+0.6*lamBarZ*nz), Cmz*(1+0.6*nz));
  }
  if(sec.isBox){ // Table B.1   not susceptible to torsional deformation
    kyz=kzz;
    kzy=Math.min(0.8*kyy, kyy);
  } else { // Table B.2   susceptible to torsional deformation (I/H/channel)
    kyz=0.6*kzz;
    const denom=Math.max(CmLT-0.25,0.05)*Math.max(chiZ,0.05);
    if(cl.cls<=2) kzy=Math.max(1-0.1*lamBarZ*nz/denom, 1-0.1*nz/denom);
    else kzy=Math.max(1-0.1*lamBarZ*nz/denom, 1-0.1*nz/denom); // same form, NA keeps the 0.1 coefficient for both class groups in this term
  }
  const McRdLT = sec.isBox? McRd : MbRd;
  const u1=ny + kyy*Mx/Math.max(McRdLT,1e-9);          // y-y (in-plane) interaction, Eq 6.61
  const u2=nz + kzy*Mx/Math.max(McRdLT,1e-9);          // z-z (out-of-plane) interaction, Eq 6.62

  const span=a.deflection?a.deflection.span:a.L, divisor=S.divisor, dlimit=span/divisor;
  const dmax=Math.abs(a.deflection?a.deflection.dmax:a.dmax), defOk=dmax<=dlimit;

  const utils=[
    {name:"Shear  Ved/Vpl,Rd",val:Fv/VplRd},
    {name:"Bending  Med/MN,Rd",val:localUtil},
    {name: sec.isBox?"Bending Med/Mc,Rd":"LTB  Med/Mb,Rd", val:ltbUtil},
    {name:"Buckling (y-y)",val:u1},
    {name:"Buckling (z-z)",val:u2},
    {name:"Deflection",val:dmax/dlimit},
  ];
  let gov=utils[0]; utils.forEach(u=>{ if(u.val>gov.val) gov=u; });
  const pass=unsupported.length===0 && utils.every(u=>u.val<=1.0001);

  return {eps,cl,clsName,unsupported,fy,Av,VplRd,Fv,lowShear,shearBuckle,McRd,hsNote,Zx,Sx,
    Ag,NRd,F,n,MNRd,awNote,Mx,localUtil,isCant,
    C1,c1method:c1r.method,LE,lamLT,Mcr,curveInfo,chiLT,kc,fmod,chiLTmod,MbRd,ltbUtil,
    lam1,curveY,curveZ,lamBarY,lamBarZ,chiY,chiZ,NbRdY,NbRdZ,Fc,Cm,cmMethod,Cmy,Cmz,CmLT,kyy,kzz,kyz,kzy,n_ratio,psiInteraction,McRdLT,u1,u2,
    span,divisor,dlimit,dmax,defOk,utils,gov,pass};
}

// ---- EN 1993-1-1 Table B.3: equivalent uniform moment factor Cm from the
// governing combination's own moment diagram (single span between restraints).
// Returns {Cm, label}. Loading type: 'uniform' / 'concentrated'; where the two
// give different values (mixed loading) the LARGER (conservative) is used.
function cmTableB3(a){
  const Mm=Math.abs(a.Mmax);
  if(Mm<1e-9) return {Cm:1.0,label:'negligible moment'};
  const fb=a.governM.fb, L=a.L;
  const M0=interpAt(fb.xs,fb.M,1e-4)/1e6, ML=interpAt(fb.xs,fb.M,L-1e-4)/1e6;
  const Ms=interpAt(fb.xs,fb.M,L/2)/1e6;
  const gfac=a.governM.combo.factors;
  const pieces=comboLoadPieces(a.governM.combo);   // pattern-aware effective load list
  let hasDist=Math.abs(gfac.G??0)>0, hasConc=false;   // auto self-weight is distributed
  pieces.forEach(p=>{ if(!p.factor) return;
    if(p.type==='point') hasConc=true; else if(p.type==='udl'||p.type==='trap') hasDist=true; });
  // 'linear end-moment diagram' means the BMD is actually a straight line
  // between the ends (no transverse-load curvature) - test every grid value
  // against the chord, not just the midpoint (a cantilever's Mmax sits
  // at the end but its diagram is far from linear; Table B.3 alpha_s applies).
  const isLinear=fb.xs.every((x,i)=>x<1e-4||x>L-1e-4||Math.abs(fb.M[i]/1e6-(M0+(ML-M0)*x/L))<=1e-5*Mm);
  if(isLinear){
    const Mh2=Math.abs(M0)>=Math.abs(ML)? M0:ML, Mo2=Math.abs(M0)>=Math.abs(ML)? ML:M0;
    const psi=Math.max(-1,Math.min(1,Mo2/(Mh2||1e-9)));
    return {Cm:Math.max(0.6+0.4*psi,0.4),label:'linear end-moment diagram, &psi; = '+psi.toFixed(2)};
  }
  // Table B.3's transverse-load diagrams do not describe arbitrary partial,
  // multiple or reversing loads. No beneficial Cm is inferred from a single
  // midpoint for those layouts.
  const active=pieces.filter(p=>Math.abs(p.factor)>1e-12);
  const simpleSpan=endsList().every(e=>e.uz)&&!(S.hinges||[]).length;   // both ends held vertically (pinned or fixed), no hinge
  const canonical=active.every(p=>p.type==='udl'&&p.x1<=1e-9&&Math.abs(p.x2-a.L)<1e-6||(p.type==='point'&&Math.abs(p.pos-a.L/2)<1e-6));
  if(!simpleSpan||!canonical||(hasDist&&hasConc)) return {Cm:1,label:'C_m = 1: arbitrary or mixed moment diagram; no Table B.3 reduction assumed'};
  const Mh=Math.abs(M0)>=Math.abs(ML)? M0:ML;
  const Mo=Math.abs(M0)>=Math.abs(ML)? ML:M0;
  const psi=Math.abs(Mh)>1e-9? Math.max(-1,Math.min(1,Mo/Mh)) : 1;
  const evalType=(type)=>{
    if(Math.abs(Mh)<0.02*Mm){ // no significant end moment: alpha_h family with Mh=0
      return type==='uniform'? 0.95 : 0.90;
    }
    if(Math.abs(Ms)<=Math.abs(Mh)+1e-12){
      const as=Math.max(-1,Math.min(1,Ms/Mh));
      if(as>=0) return 0.2+0.8*as;
      if(psi>=0) return type==='uniform'? 0.1-0.8*as : -0.8*as;
      return type==='uniform'? 0.1*(1-psi)-0.8*as : 0.2*(-psi)-0.8*as;
    }
    const ah=Math.max(-1,Math.min(1,Mh/Ms));
    const mod=(ah<0&&psi<0)? (1+2*psi) : 1;
    return type==='uniform'? 0.95+0.05*ah*mod : 0.90+0.10*ah*mod;
  };
  let Cm,label;
  const tag=Math.abs(Ms)<=Math.abs(Mh)? '&alpha;<sub>s</sub> = '+ (Math.abs(Mh)<0.02*Mm?'&mdash;':(Ms/Mh).toFixed(3)) : '&alpha;<sub>h</sub> = '+(Mh/Ms).toFixed(3);
  if(hasDist&&hasConc){ Cm=Math.max(evalType('uniform'),evalType('concentrated')); label='mixed loading (larger of uniform/concentrated), '+tag; }
  else if(hasConc){ Cm=evalType('concentrated'); label='concentrated load diagram, '+tag; }
  else { Cm=evalType('uniform'); label='uniform load diagram, '+tag; }
  return {Cm:Math.max(Cm,0.4),label:label+', M<sub>h</sub> = '+Mh.toFixed(1)+', M<sub>s</sub> = '+Ms.toFixed(1)+', &psi; = '+psi.toFixed(2)};
}
// ---- EN 1993-1-1 cl 6.3.3, Annex B Method 2 (interaction factors kyy, kzz,
// kyz, kzy from Table B.1 [not susceptible to torsional deformation] or B.2
// [susceptible]; Cm from Table B.3). This is a compression-member check, so
// callers should only create it when N_Ed is compressive.
function lcrZFromRestraints(a){
  // SCI P360 Section 6.2: secondary members (purlins, side rails, ties) attached
  // to a member "act as bracing points" and reduce the MINOR-AXIS (z-z) strut
  // buckling length; between adjacent lateral restraints k = 1.0 and L is the
  // distance between restraint points. Applied only in the EC3 unrestrained
  // mode, where the intermediate-lateral-restraint list is visible and edited;
  // only restraints that hold lateral displacement v count; the bays run
  // between the ends that hold U_y. Returns null (keep the end-fixity length)
  // when there are no intermediate v-restraints, when only one end holds U_y
  // (lateral cantilever), or if anything is degenerate.
  if(!(S.code==='EC3' && (S.restraint||'full')!=='full')) return null;
  const ir=(S.ltbRestraints||[]).filter(r=>r.v!==false).map(r=>(+r.pos)*1000).filter(x=>isFinite(x)&&x>=0&&x<=a.L);
  if(!ir.length) return null;
  const endPts=endsList().filter(e=>e.uy).map(e=>e.x*1000);   // ends holding lateral translation U_y
  if(endPts.length<2) return null;               // lateral cantilever: keep the end-fixity default
  const pts=[...new Set(endPts.concat(ir).map(x=>+x.toFixed(3)))].sort((p,q)=>p-q);
  if(pts.length<2) return null;
  let gmax=0;
  for(let i=1;i<pts.length;i++) gmax=Math.max(gmax,pts[i]-pts[i-1]);
  return gmax>1e-6? gmax : null;
}
/* Torsional buckling length from the TWIST restraints (19 Sep 2026 review
   finding): EN 1993-1-1 6.3.1.4(1) with EN 1993-1-3 6.2.3(5) - l_T is set by
   the torsional / warping restraint at the ends of the torsional segment. A
   restraint that holds lateral displacement v only (phi unticked) does not
   bound the torsional mode, so L_T is the largest spacing between points that
   prevent twist: every end whose R_x is restrained and every intermediate
   restraint with phi !== false. Mirrors lcrZFromRestraints(): null (keep the
   flexural default) when there are no intermediate twist restraints, when
   only one end holds twist (torsion cantilever: L_cr,y), or in the fully
   restrained mode. */
function lcrTFromTwistRestraints(a){
  if(!(S.code==='EC3' && (S.restraint||'full')!=='full')) return null;
  const ir=(S.ltbRestraints||[]).filter(r=>r.phi!==false).map(r=>(+r.pos)*1000).filter(x=>isFinite(x)&&x>=0&&x<=a.L);
  if(!ir.length) return null;
  const endPts=endsList().filter(e=>e.rx).map(e=>e.x*1000);
  if(endPts.length<2) return null;
  const pts=[...new Set(endPts.concat(ir).map(x=>+x.toFixed(3)))].sort((p,q)=>p-q);
  if(pts.length<2) return null;
  let gmax=0;
  for(let i=1;i<pts.length;i++) gmax=Math.max(gmax,pts[i]-pts[i-1]);
  return gmax>1e-6? gmax : null;
}
/* ---- k_c floor (19 Sep 2026 gap closure, item 3.5) ----
   NA 2.18 allows k_c = 1/sqrt(C1). Table 6.6 lists k_c down to 0.60 (the
   psi = -1 end-moment case, C1 = 2.76), so a back-calculated or Serna C1 above
   2.76 is capped at that lower bound: k_c >= 1/sqrt(2.76) = 0.602. */
const KC_FLOOR_C1=2.76, KC_FLOOR=1/Math.sqrt(KC_FLOOR_C1);
function kcFromC1(C1){
  const raw=Math.min(1/Math.sqrt(Math.max(C1,1e-6)),1);
  return {kc:Math.max(raw,KC_FLOOR), kcRaw:raw, floored:raw<KC_FLOOR-1e-12};
}
/* ---- Effective area of a Class-4 web in uniform compression (item 1.9(c) /
   3.9(c)), EN 1993-1-5 4.4: internal compression element, psi = 1, k_sigma = 4,
   lambda_p = (b/t)/(28.4 eps sqrt(k_sigma)), rho = (lambda_p - 0.055(3 + psi))/
   lambda_p^2 <= 1 (lambda_p > 0.673; = 1 otherwise), b_eff = rho b split
   b_e1 = b_e2 = 0.5 b_eff. b = the flat web depth of Table 5.2 (sec.d: I/H
   h - 2t_f - 2r, hollow h - 3t); a hollow section has two such walls. The
   reduction is symmetric, so the centroid does not move (e_N = 0). Applies
   when the web is Class 4 in uniform compression (d/t > 42 eps); the flanges
   must stay Class <= 3 in compression (outstand 14 eps / internal 42 eps) -
   the caller blocks otherwise. Pure. */
function aeffWebCompression(sec,eps){
  const A=sec.A*1e2, tw=sec.tw, bbar=sec.d, ratio=sec.dt, nWebs=sec.isBox? 2 : 1;
  const out={applies:false,A,Aeff:A,ratio:1,nWebs,bbar,tw,eps,limit:42*eps,dt:ratio,ksig:4,psi:1,lamP:null,rho:1,beff:bbar,be1:bbar/2,be2:bbar/2,bineff:0,eN:0};
  if(!(ratio>42*eps) || !(bbar>0&&tw>0)) return out;
  const lamP=ratio/(28.4*eps*Math.sqrt(4));
  const rho= lamP<=0.673? 1 : Math.min((lamP-0.055*(3+1))/(lamP*lamP),1);
  const beff=rho*bbar, bineff=(1-rho)*bbar;
  const Aeff=A-nWebs*bineff*tw;
  return Object.assign(out,{applies:true,lamP,rho,beff,be1:beff/2,be2:beff/2,bineff,Aeff,ratio:Aeff/A});
}
/* ---- Flange outstand stresses of an I/H section under N + M_y + M_z (item
   1.9(b), printed with the classification): elastic extreme-fibre values,
   compression positive. The outstand of the compression flange on the side
   compressed by M_z has root stress sN + sMy + sMz(root) and tip stress sN +
   sMy + sMz(tip): both compressive -> alpha = 1 -> the uniform-compression
   bound 9e/10e/14e governs. The opposite outstand is relieved at its tip
   (tip in tension when sMz(tip) > sN + sMy) and takes the laxer stress-
   gradient limits of Table 5.2 sheet 2, so it never governs. Pure. */
function mzFlangeStress(sec,F,My,Mz){
  const A=sec.A*1e2, Zx=sec.Zx*1e3, Zy=sec.Zy*1e3;
  const sN=Math.max(F,0)*1000/A, sMy=Math.abs(My)*1e6/Zx, sMzTip=Math.abs(Mz)*1e6/Zy;
  const yRoot=sec.tw/2+sec.r, yTip=sec.B/2;
  const sMzRoot=sMzTip*yRoot/yTip;
  const compRoot=sN+sMy+sMzRoot, compTip=sN+sMy+sMzTip;
  const relRoot=sN+sMy-sMzRoot, relTip=sN+sMy-sMzTip;
  return {sN,sMy,sMzTip,sMzRoot,yRoot,yTip,compRoot,compTip,relRoot,relTip,
    compAlpha:(compRoot>=0&&compTip>=0)? 1 : null, relState: relTip<0? 'tip in tension' : 'wholly in compression'};
}
/* ---- Torsional and torsional-flexural buckling of a channel under N (item
   3.10, EN 1993-1-1 6.3.1.4). Monosymmetric about y-y (the major axis, the
   axis of symmetry); the shear centre lies on it at y0 = e_sc from the
   centroid (SCI P385 Table A.3 / section tables, sec.tp.esc, mm). i0^2 = i_y^2 +
   i_z^2 + y0^2; N_cr,T = (G I_T + pi^2 E I_w/L_T^2)/i0^2 with L_T = the
   spacing of the TWIST restraints (lcrTFromTwistRestraints: supports and
   intermediate restraints with phi ticked, never the v-only spacing L_cr,z),
   capped at L_cr,y, unless the user enters S.LT (m); the torsional mode
   couples with flexure
   about the axis of symmetry (the y-y flexural mode, N_cr,y) through the
   standard cubic, which for one axis of symmetry reduces to
   N_cr,TF = (N_cr,y + N_cr,T)/(2 beta) [1 - sqrt(1 - 4 beta N_cr,y N_cr,T/
   (N_cr,y + N_cr,T)^2)], beta = 1 - (y0/i0)^2 (EN 1993-1-3 6.2.3(6) form).
   N_cr = min(N_cr,T, N_cr,TF); lambda_T = sqrt(A f_y/N_cr); chi from the curve
   related to the z-z axis (6.3.1.4(3)): Table 6.2 "U-sections: any axis ->
   curve c" [verify]; N_b,T,Rd = chi A f_y/gamma_M1. Pure. */
function torsionalFlexuralBuckling(sec,fy,E,LcrY,LcrZ,Ag,gM1,LcrT){
  const G=81000;
  const IT=((sec.tp&&sec.tp.IT)? sec.tp.IT : sec.J)*1e4, ITSrc=(sec.tp&&sec.tp.IT)? 'SCI P385 Table A.3' : 'section table';
  const Iw=(((sec.tp&&sec.tp.Iw!=null)? sec.tp.Iw : sec.Iw)||0)*1e12;
  const y0=(sec.tp&&sec.tp.esc!=null&&isFinite(+sec.tp.esc))? +sec.tp.esc : null;
  if(y0==null || !(IT>0)) return {ok:false,reason:'the shear-centre offset e<sub>sc</sub> (SCI P385 Table A.3) or I<sub>T</sub> is not tabulated for this section, so N<sub>cr,T</sub> / N<sub>cr,TF</sub> cannot be formed'};
  const iy=sec.rx*10, iz=sec.ry*10;
  const i0sq=iy*iy+iz*iz+y0*y0, i0=Math.sqrt(i0sq);
  const LTin=(S.LT!=null && S.LT!=='' && isFinite(+S.LT) && +S.LT>0)? (+S.LT)*1000 : null;
  // default L_T: the twist-restraint spacing (LcrT, capped at L_cr,y by the
  // caller), never the lateral-only spacing L_cr,z; LcrZ is kept for reporting
  const LTdef= (LcrT!=null && isFinite(LcrT) && LcrT>0)? LcrT : LcrY;
  const LT= LTin!=null? LTin : LTdef;
  const LTSrc= LTin!=null? 'user L<sub>T</sub>' : (LcrT!=null && LcrT<LcrY-1e-6)? 'spacing of twist restraints (ends with R<sub>x</sub> held and restraints with &phi; held)' : 'L<sub>cr,y</sub> (no intermediate twist restraint)';
  const NcrT=(G*IT+Math.PI*Math.PI*E*Iw/(LT*LT))/i0sq;            // N
  const NcrY=Math.PI*Math.PI*E*sec.Ix*1e4/(LcrY*LcrY);           // N, flexural about the axis of symmetry (y-y, major)
  const beta=1-(y0*y0)/i0sq;
  const NcrTF=(NcrY+NcrT)/(2*beta)*(1-Math.sqrt(Math.max(1-4*beta*NcrY*NcrT/Math.pow(NcrY+NcrT,2),0)));
  const Ncr=Math.min(NcrT,NcrTF), mode= NcrTF<NcrT? 'TF' : 'T';
  const lamT=Math.sqrt(Ag*fy/Ncr);
  const cvT=strutCurveEC3(sec,'z');
  const chiT=chiStrutEC3(lamT,cvT.alpha);
  const NbT=chiT*Ag*fy/gM1/1000;                                 // kN
  return {ok:true,y0,y0Src:'e<sub>sc</sub> = shear centre to centroid, SCI P385 Table A.3 (section tables e<sub>0</sub> + c<sub>y</sub>)',iy,iz,i0,i0sq,IT,ITSrc,Iw,G,LT,LTSrc,LTdef,LcrZ,LcrT,
    NcrT:NcrT/1000,NcrY:NcrY/1000,beta,NcrTF:NcrTF/1000,Ncr:Ncr/1000,mode,lamT,cvT,chiT,NbT};
}
function annexB2(a,sec,fy,cl,MbRdI,useB1,isCant,aeff){
  if(a.ulsResults&&a.ulsResults.length>1){
    const rows=a.ulsResults.map(res=>({...annexB2(analysisForCombination(a,res),sec,fy,cl,MbRdI,useB1,isCant,aeff),combo:res.combo.label}));
    return rows.reduce((p,r)=>Math.max(r.u1,r.u2)>Math.max(p.u1,p.u2)?r:p);
  }
  const gM1=1.0, E=a.E;
  const Fc=Math.max(S.axial||0,0);
  const Ag=sec.A*1e2;
  // Class-4 web in uniform compression (item 1.9(c)): A_eff replaces A in
  // N_Rk, in lambda-bar (6.3.1.3(1): lambda = sqrt(A_eff f_y/N_cr)) and in the
  // Table 6.7 Class-4 column of Eq 6.61/6.62 (W_eff,y = W_el,y, e_N = 0).
  const aeffOn=!!(aeff&&aeff.applies&&Fc>1e-9);
  const Aeff= aeffOn? aeff.Aeff : Ag;
  const aeffFac= aeffOn? Math.sqrt(Aeff/Ag) : 1;
  // Strut lengths per axis (destabilising x1.2 is an LTB concept, not applied).
  // Both axes default from the END FIXITIES (lcrDefaults, js/03-state-ui.js:
  // P360 Table 6.2 / BS 5950 Table 22 style - 0.7 L fixed-fixed, 0.85 L one
  // end fixed, 1.0 L pinned-pinned, 2.0 L fixed-free (cantilever), 1.2 L
  // fixed-guided; y-y from U_z / R_y, z-z from U_y / R_z, printed [verify]);
  // the entered L_E/L factor overrides both. Minor axis: further reduced to
  // the largest spacing between adjacent lateral restraint points where
  // intermediate restraints are modelled (SCI P360 6.2 bracing-point
  // assumption), never longer than its end-fixity length.
  const lcr=lcrDefaults(S);
  if(Fc>1e-9 && (lcr.Ky==null || lcr.Kz==null)) throw 'Strut buckling: the end fixities form a mechanism ('+(lcr.Ky==null? lcr.basisY : lcr.basisZ)+'); N_Ed cannot be carried.';
  const cantStrut = !!isCant && Fc>1e-9;
  const leOverride = !!lcr.override;
  const Ky = lcr.Ky!=null? lcr.Ky : 1.0;
  const KzEnd = lcr.Kz!=null? lcr.Kz : 1.0;
  const lcrBasis = lcr.basis;
  const LcrY=Ky*a.L;
  const LcrZEnd=KzEnd*a.L;
  const lz=lcrZFromRestraints(a);
  const LcrZ=lz!=null? Math.min(lz,LcrZEnd) : LcrZEnd;
  const Kz=LcrZ/a.L;
  const lczFromRestraints=(lz!=null && LcrZ<LcrZEnd-1e-6);
  const Lcr=LcrY;                                 // kept for report compatibility
  const lam1=Math.PI*Math.sqrt(E/fy);
  const rx=sec.rx*10, ry=sec.ry*10;
  const lamY=(LcrY/rx)/lam1*aeffFac, lamZ=(LcrZ/ry)/lam1*aeffFac;
  const cvY=strutCurveEC3(sec,'y'), cvZ=strutCurveEC3(sec,'z');
  const chiY=chiStrutEC3(lamY,cvY.alpha), chiZ=chiStrutEC3(lamZ,cvZ.alpha);
  const NbY=chiY*Aeff*fy/gM1/1000, NbZ=chiZ*Aeff*fy/gM1/1000;   // kN
  // channel under compression: torsional / torsional-flexural buckling (6.3.1.4);
  // the lower of chi_T and the flexural chi feeds both axial terms of 6.61/6.62
  let tfb=null;
  if(sec.kind==='channel' && Fc>1e-9){
    // L_T from the twist-restraint spacing (phi held), capped at L_cr,y; a
    // lateral-only restraint shortens L_cr,z but not L_T
    const lt=lcrTFromTwistRestraints(a);
    const LcrT= lt!=null? Math.min(lt,LcrY) : LcrY;
    tfb=torsionalFlexuralBuckling(sec,fy,E,LcrY,LcrZ,Aeff,gM1,LcrT);
    if(tfb.ok) tfb.util=Fc/Math.max(tfb.NbT,1e-9);
  }
  const NbYeff= (tfb&&tfb.ok)? Math.min(NbY,tfb.NbT) : NbY;
  const NbZeff= (tfb&&tfb.ok)? Math.min(NbZ,tfb.NbT) : NbZ;
  const ny=Fc/Math.max(NbYeff,1e-9), nz=Fc/Math.max(NbZeff,1e-9);
  const cm=cmTableB3(a);
  let Cmy=cm.Cm, CmLT=cm.Cm, swayNote=false;
  if(isCant && Fc>1e-6 && Cmy<0.9){ Cmy=0.9; CmLT=0.9; swayNote=true; } // Table B.3 note: sway buckling mode -> Cm = 0.9
  const Cmz=1.0;                                   // 20 Sep 2026 review: Table B.3 upper bound - exact for the imposed constant M_z (psi = 1), conservative for the twist-induced phi.M_y diagram (M_z,Ed = M_z + phi.M_y below)
  // k_ij column: Class 1/2 plastic forms only for a Class 1/2 section that is not
  // a channel (elastic column, conservative) and has no Class-4 web in
  // compression (Table 6.7 Class-4 column -> the Class 3/4 rows of Table B.1/B.2)
  const c12=cl.cls<=2 && sec.kind!=='channel' && !aeffOn;
  const rhsRow=!!sec.isBox;                        // Table B.1 "rectangular hollow sections" row for k_zz (item 3.12(a))
  const kyy = c12? Math.min(Cmy*(1+(lamY-0.2)*ny), Cmy*(1+0.8*ny))
                 : Math.min(Cmy*(1+0.6*lamY*ny),   Cmy*(1+0.6*ny));
  const kzz = c12? (rhsRow? Math.min(Cmz*(1+(lamZ-0.2)*nz), Cmz*(1+0.8*nz))       // RHS row: same form as k_yy
                          : Math.min(Cmz*(1+(2*lamZ-0.6)*nz), Cmz*(1+1.4*nz)))    // I-section row
                 : Math.min(Cmz*(1+0.6*lamZ*nz),     Cmz*(1+0.6*nz));
  const kyz = (c12&&!rhsRow)? 0.6*kzz : kzz;      // Table B.1: I-sections Class 1/2 0.6k_zz; RHS and the Class 3 column k_zz (shared by B.2)
  let kzy, kzyLbl;
  if(useB1){ kzy=(c12?0.6:0.8)*kyy; kzyLbl='Table B.1: '+(c12?'0.6':'0.8')+'k<sub>yy</sub>'; }
  else {
    const dn=Math.max(CmLT-0.25,1e-6), lz=Math.min(lamZ,1);
    const coef=c12? 0.1 : 0.05;
    let e1=1-coef*lz*nz/dn;
    if(c12 && lamZ<0.4) e1=Math.min(0.6+lamZ,e1);
    kzy=Math.max(e1,0); kzyLbl='Table B.2 (susceptible): 1&minus;'+coef+'&middot;min(&lambda;&#772;<sub>z</sub>,1)&middot;n<sub>z</sub>/(C<sub>mLT</sub>&minus;0.25)';
  }
  // Table 6.7 Class-4 column (A_eff case): M_y,Rk = W_eff,y f_y with W_eff,y =
  // W_el,y (flanges Class <= 3 in compression, web Class 4 only in uniform
  // compression, so the section under bending keeps its elastic modulus), so a
  // Class 1/2 M_b,Rd handed in (plastic W_y) is scaled by W_el,y/W_pl,y.
  const wFac= (aeffOn && cl.cls<=2)? sec.Zx/sec.Sx : 1;
  const Mrd=Math.max(MbRdI*wFac,1e-9), Mx=Math.abs(a.Mmax);
  // Minor-axis moment terms (fully biaxial 6.61/6.62). Mz,Ed is the direct design
  // input; its resistance Mc,z,Rd carries no LTB reduction (chi_LT is major-axis only).
  // 20 Sep 2026 torsion + N/Mz: M_z,Ed = imposed M_z + max over the member of |phi.M_y|,
  // the second-order minor-axis moment of this combination's torsion solution (EN 1993-1-1
  // 5.2.1(3); SCI P385 3.1.2 M_z = phi.M_y); zero without torsion, so nothing changes then.
  const MzImp=Math.abs(S.Mz||0);
  const cb0=(a.ulsResults&&a.ulsResults[0])? a.ulsResults[0].combo : (a.governM? a.governM.combo : null);
  const MzTwist=cb0? torsionMzTwistMax(a,cb0) : 0;
  const MzEd=MzImp+MzTwist;
  const Mcz=Math.max((c12? sec.Sy : sec.Zy)*1e3*fy/gM1/1e6, 1e-9);   // Mc,z,Rd
  const mzTerm = MzEd>1e-9 ? MzEd/Mcz : 0;                            // Mz,Ed / Mc,z,Rd
  const u1=ny + kyy*Mx/Mrd + kyz*mzTerm;           // Eq 6.61
  const u2=nz + kzy*Mx/Mrd + kzz*mzTerm;           // Eq 6.62
  return {Fc,Mx,Lcr,LcrY,LcrZ,Ky,Kz,KzEnd,cantStrut,leOverride,lcrBasis,lcrBasisY:lcr.basisY,lcrBasisZ:lcr.basisZ,lczFromRestraints,lam1,lamY,lamZ,cvY,cvZ,chiY,chiZ,NbY,NbZ,NbYeff,NbZeff,ny,nz,
    Cmy,Cmz,CmLT,cmLabel:cm.label,swayNote,useB1,c12,rhsRow,kyy,kzz,kyz,kzy,kzyLbl,MbRdI,MbRdEff:Mrd,wFac,Mcz,MzEd,MzImp,MzTwist,mzTerm,biax:MzImp>1e-9,u1,u2,
    aeffOn,Aeff,Ag,aeffFac,tfb};
}
/* ---------------------------------------------------------------------------
   EN 1993-1-5 clause 6: resistance of the web to transverse forces, with the
   clause 7.2 interaction (19 Sep 2026 gap closure, items 2.16 + 2.18). Pure:
   reads the analysis a (every ULS combination's own reactions and diagram),
   the section, fy, eps, the class and the state S (per-load ss / stiff, per-
   support ss / stiff, axial). No DOM.
   Stations: every point load (non-zero P) and every support. Loads within
   0.5 mm of each other share a station (forces summed, the smaller s_s
   kept); a point load at a support position is a "both" station.
   Load types (Figure 6.1): (a) interior load resisted by shear in the web,
   k_F = 6 + 2(h_w/a)^2; (b) load transferred through the web to the opposite
   flange (a point load directly over a support), k_F = 3.5 + 2(h_w/a)^2,
   F_Ed = max(P, R) taken conservatively as the through-load; (c) load near an
   unstiffened end, k_F = 2 + 6(s_s + c)/h_w <= 6, evaluated whenever
   (s_s + c) < 2h_w/3 (the value at which k_F(c) reaches the long-panel 6),
   together with type (a); the lower F_Rd governs.
   F_Rd = f_yw L_eff t_w / gamma_M1 (6.2), L_eff = chi_F l_y, chi_F = 0.5 /
   lambda_F <= 1, lambda_F = sqrt(l_y t_w f_yw / F_cr), F_cr = 0.9 k_F E t_w^3
   / h_w (6.4); m1 = f_yf b_f/(f_yw t_w), m2 = 0.02 (h_w/t_f)^2 if lambda_F >
   0.5 else 0 (one re-evaluation with m2 = 0 when the first pass gives
   lambda_F <= 0.5); l_y = s_s + 2 t_f (1 + sqrt(m1 + m2)) <= a for (a)/(b),
   the two 6.5(4) expressions with l_e = k_F E t_w^2/(2 f_yw h_w) <= s_s + c
   for (c). b_f is limited to 15 eps t_f each side of the web (6.5(1)).
   s_s: per-load input (default 0) and per-support input (blank = the LOWER
   BOUND s_s = 0, since the seating length is set by the bearing, not by the
   beam; F_Rd rises monotonically with s_s, so a station that passes at 0 is
   verified for any seating, while a station that fails at 0 with s_s not
   entered is reported NOT VERIFIED - blocking, with the s_s = 0 values
   printed - instead of FAIL: 19 Sep 2026 review finding, replacing the
   former default s_s = B), capped at h_w (6.3(1)). c = max(d - s_s/2, 0) with d the distance
   from the station to the nearer member end. a = distance between declared
   bearing stiffeners bounding the station, the full member length when none
   are declared (conservative, printed).
   Hollow sections: two webs, each a plate of thickness t with the tabulated
   flat depth (h_w = d, corner geometry of the section table); the flange
   width per web is B/2 limited to t + 15 eps t; the load is shared between
   the webs by the lever rule of its eccentricity e (0.5 each when e = 0).
   Interaction 7.2: eta_2 + 0.8 eta_1 <= 1.4 with eta_2 = F_Ed/F_Rd and eta_1
   = M_Ed/M_c,Rd (+ N_Ed/N_pl,Rd when N_Ed != 0) at the same station in the
   same combination; M_c,Rd is the unreduced class-consistent W_y f_y/gamma_M0
   [verify: EN 1993-1-5 4.6 writes eta_1 with W_eff, i.e. W_el for Class 3].
   The interaction is evaluated at every station; where the loaded flange is
   in tension 7.2(2) refers to 6.2.1(5) of EN 1993-1-1 instead and the row
   says so (the 7.2 expression is still applied as a screen).
   A station whose "bearing stiffener provided" box is ticked is not checked:
   it prints "stiffener declared - design stiffener separately (EN 1993-1-5
   9.4)" as an advisory. Returns {checked, stations, util2, util72, ...}.
   --------------------------------------------------------------------------- */
const WEB_UTIL_NAMES=["Web transverse force  F_Ed/F_Rd (EN 1993-1-5 6.2)","Web transverse force + bending (EN 1993-1-5 7.2)"];
const TFB_UTIL_NAME="Torsional-flexural buckling  N_Ed/N_b,T,Rd (6.3.1.4)";
function mvnUtilName(m){ return "Bending+shear+"+(m.N>1e-9? (m.biax? "axial+biaxial" : "axial") : "biaxial")+" (6.2.10)"; }
function webTransverseCheck(a,sec,fy,eps,cl){
  const gM0=1.0, gM1=1.0, E=a.E, L=a.L;
  const isBox=!!sec.isBox, chan=sec.kind==='channel';
  const tw=sec.tw, tf=sec.tf;
  const hw= isBox? sec.d : sec.D-2*sec.tf;
  const nWebs= isBox? 2 : 1;
  const fyw=fy, fyf=fy;
  const bfRaw= isBox? sec.B/2 : sec.B;
  const bfLim= (isBox||chan)? tw+15*eps*tf : tw+30*eps*tf;
  const bf=Math.min(bfRaw,bfLim);
  const m1=fyf*bf/(fyw*tw);
  const m2full=0.02*Math.pow(hw/tf,2);
  const Wy=(cl.cls<=2? sec.Sx : sec.Zx)*1e3;
  const McRd0=Wy*fy/gM0/1e6;                     // kN.m, unreduced, class-consistent
  const NplRd=sec.A*100*fy/gM0/1000;             // kN
  const NEd=Math.abs(S.axial||0);
  const out={checked:false,hw,tw,tf,bf,bfRaw,bfLim,m1,m2full,eps,fyw,fyf,nWebs,isBox,chan,McRd0,NplRd,NEd,cls:cl.cls,
    stations:[],gov2:null,gov72:null,util2:0,util72:0,unsupported:[],advisory:[],aBasis:'',anyDefaultSs:false,anyTension:false};
  if(!(hw>0&&tw>0&&tf>0)){ out.unsupported.push('Web transverse forces (EN 1993-1-5 clause 6): the web geometry (h_w, t_w, t_f) is undefined for this section; the check cannot be made.'); return out; }
  const tol=0.5;   // mm: loads / supports closer than this share a station
  const st=[];
  const find=x=>st.find(s=>Math.abs(s.x-x)<=tol);
  const ssOf=v=>(v!=null && v!=='' && Number.isFinite(+v))? Math.max(+v,0) : null;
  S.loads.forEach((ld,i)=>{
    if(ld.isSelfWeight||ld.type!=='point'||!(Math.abs(+ld.P||0)>1e-12)) return;
    const x=(+ld.pos)*1000;
    let s=find(x); if(!s){ s={x,loads:[],support:null}; st.push(s); }
    s.loads.push({i,ss:ssOf(ld.ss),stiff:!!ld.stiff,e:(S.eccOn&&Number.isFinite(+ld.e))? Math.abs(+ld.e) : 0});
  });
  // every end that carries a vertical reaction (U_z held): i = index into the
  // solver's reaction list, n = end number
  verticalEnds(S).forEach(sp=>{
    const x=(+sp.pos)*1000;
    let s=find(x); if(!s){ s={x,loads:[],support:null}; st.push(s); }
    s.support={i:sp.i,n:sp.end,ss:ssOf(sp.ss),stiff:!!sp.stiff,type:sp.type};
  });
  st.sort((p,q)=>p.x-q.x);
  // declared stiffeners bound the web panels (a); none declared -> a = L
  const stiffX=st.filter(s=>s.loads.some(l=>l.stiff)||(s.support&&s.support.stiff)).map(s=>s.x);
  out.aBasis= stiffX.length? 'a = distance between the declared bearing stiffeners bounding the station (member ends otherwise)' : 'no transverse stiffeners declared: a = L = '+g(L,0)+' mm, the full member length (conservative)';
  const panelOf=x=>{ let l=0,r=L; stiffX.forEach(p=>{ if(p<x-tol && p>l) l=p; if(p>x+tol && p<r) r=p; }); return {a:Math.max(r-l,1e-6),l,r}; };
  const solve=(type,ss,c,aPanel)=>{
    const kF= type==='a'? 6+2*Math.pow(hw/aPanel,2) : type==='b'? 3.5+2*Math.pow(hw/aPanel,2) : Math.min(2+6*(ss+c)/hw,6);
    const Fcr=0.9*kF*E*Math.pow(tw,3)/hw;                       // N
    const leRaw= type==='c'? kF*E*tw*tw/(2*fyw*hw) : null;
    const le= type==='c'? Math.min(leRaw,ss+c) : null;
    const lyFor=(m2)=>{
      if(type==='c'){ const l1=le+tf*Math.sqrt(m1/2+Math.pow(le/tf,2)+m2), l2=le+tf*Math.sqrt(m1+m2); return {ly:Math.min(l1,l2),l1,l2,capA:false}; }
      const l=ss+2*tf*(1+Math.sqrt(m1+m2)); return {ly:Math.min(l,aPanel),l1:l,l2:null,capA:l>aPanel};
    };
    let m2=m2full, r=lyFor(m2), lam=Math.sqrt(r.ly*tw*fyw/Fcr), iter=false, lam1=lam, ly1=r.ly;
    if(lam<=0.5){ m2=0; r=lyFor(0); lam=Math.sqrt(r.ly*tw*fyw/Fcr); iter=true; }
    const chiRaw=0.5/lam, chi=Math.min(chiRaw,1);
    const Leff=chi*r.ly;
    const FRd=fyw*Leff*tw/gM1/1000;                             // kN per web
    return {type,kF,Fcr:Fcr/1000,leRaw,le,m2,iter,lam1,ly1,ly:r.ly,l1:r.l1,l2:r.l2,capA:r.capA,lam,chiRaw,chi,Leff,FRd};
  };
  st.forEach(s=>{
    const kind= s.loads.length&&s.support? 'both' : s.support? 'support' : 'load';
    const stiff= s.loads.some(l=>l.stiff)||(s.support&&s.support.stiff);
    const n= s.support? s.support.n : null;
    const loadIdx=s.loads.map(l=>l.i+1);
    const label= kind==='support'? 'End '+n+' reaction' : kind==='load'? 'point load '+loadIdx.join('+') : 'point load '+loadIdx.join('+')+' over End '+n;
    const rec={x:s.x,kind,label,n,loadIdx,stiff:!!stiff,isEnd:(s.x<=tol||s.x>=L-tol)};
    if(stiff){
      rec.msg='stiffener declared - design stiffener separately (EN 1993-1-5 9.4)';
      out.advisory.push('Web transverse forces at x = '+g(s.x/1000,3)+' m ('+label+'): '+rec.msg+'.');
      out.stations.push(rec); return;
    }
    // stiff bearing length: the smaller of the entries at the station; a blank
    // support entry is the lower bound 0 (ssDefault: NOT VERIFIED if it fails)
    const ssLoad= s.loads.length? Math.min(...s.loads.map(l=>l.ss==null? 0 : l.ss)) : null;
    let ssSup=null, ssDefault=false;
    if(s.support){ ssSup= s.support.ss==null? 0 : s.support.ss; ssDefault= s.support.ss==null; }
    const ssIn= kind==='both'? Math.min(ssLoad,ssSup) : kind==='load'? ssLoad : ssSup;
    const ssCap= ssIn>hw;
    const ss=Math.min(ssIn,hw);
    if(ssDefault) out.anyDefaultSs=true;
    const d=Math.min(s.x,L-s.x);
    const c=Math.max(d-ss/2,0);
    const endZone=(ss+c)<2*hw/3;
    const panel=panelOf(s.x);
    const types= kind==='both'? ['b'].concat(endZone? ['c']:[]) : (endZone? ['a','c'] : ['a']);
    const sols=types.map(t=>solve(t,ss,c,panel.a));
    const gov=sols.reduce((p,q)=>q.FRd<p.FRd? q : p);
    // load share per web (box: lever rule of the largest eccentricity at the station)
    const eMax=s.loads.length? Math.max(...s.loads.map(l=>l.e)) : 0;
    const share= isBox? Math.min(1,0.5+eMax/Math.max(sec.B-tw,1e-9)) : 1;
    const FRdTot= isBox? gov.FRd/share : gov.FRd;
    // F_Ed, M_Ed per ULS combination (its own load pieces, reactions and diagram);
    // the gamma_G,inf = 1.0 STR set-B companions (19 Sep 2026 review) are swept
    // too - a relieving G raises the reaction where G lifts the support
    const ulsList=a.ulsResults.concat((a.ulsCompanions||[]).filter(res=>res.combo.gInfSet==='B'));
    const cases=ulsList.map(res=>{
      let P=0;
      comboLoadPieces(res.combo).forEach(p=>{ if(p.type==='point' && Math.abs(p.pos-s.x)<=tol) P+=p.P*p.factor; });
      const R= s.support? Math.max(res.r.reactions[s.support.i].V/1000,0) : 0;
      const F= kind==='both'? Math.max(Math.abs(P),R) : kind==='load'? Math.abs(P) : R;
      // M_Ed at the station: an end station reads the diagram a fraction
      // inside the member (x = 1e-4 mm, a grid point of sfdBmd, as analyse()
      // does for M0end / MLend) - the reaction couple of a fixed / guided end
      // is inside it, whereas sfdBmd closes the diagram to zero at x = L
      // beyond the End 2 reaction couple, so a sample exactly at the station
      // dropped the hogging end moment from the 7.2 interaction (19 Sep 2026
      // review finding F-A); an interior station reads the larger-magnitude
      // side of a jump (mAtStation: an applied couple at the same x)
      const xEnd= s.x<=tol? 1e-4 : s.x>=L-tol? L-1e-4 : null;
      const Ms=(xEnd!=null? interpAt(res.fb.xs,res.fb.M,xEnd) : mAtStation(res.fb,s.x,0,L))/1e6;
      const M=Math.abs(Ms);
      const flange= kind==='support'? 'bottom' : (P>=0? 'top' : 'bottom');
      // is the loaded flange the compression flange (7.2(1)) - a "both" station loads both flanges,
      // and a station without a coincident moment (simply supported end: the sample a fraction
      // inside the member is R x 1e-4 mm, negligible against the diagram peak) has no tension flange
      const noM = M<1e-5*Math.max(Math.abs(res.Mmax)/1e6,1e-9);
      const flangeComp= kind==='both'? true : noM? true : (flange==='top'? Ms>1e-9 : Ms<-1e-9);
      const flangeState= kind==='both'? 'load through the web' : noM? 'no coincident moment' : (flangeComp? 'in compression' : 'in tension');
      const eta2=F*share/Math.max(gov.FRd,1e-9);
      const eta1=M/Math.max(McRd0,1e-9)+(NEd>1e-9? NEd/Math.max(NplRd,1e-9) : 0);
      const u72raw=eta2+0.8*eta1, u72=u72raw/1.4;
      return {combo:res.combo.label,P,R,F,Fweb:F*share,M,Ms,flange,flangeComp,flangeState,eta2,eta1,u72raw,u72};
    });
    let g2=0,g72=0; cases.forEach((cs,i)=>{ if(cs.eta2>cases[g2].eta2) g2=i; if(cs.u72>cases[g72].u72) g72=i; });
    if(!cases[g2].flangeComp && cases[g2].F>1e-9) out.anyTension=true;
    Object.assign(rec,{ssIn,ss,ssCap,ssDefault,d,c,endZone,a:panel.a,panel,types,sols,gov,type:gov.type,share,eMax,FRd:gov.FRd,FRdTot,cases,g2,g72,
      eta2:cases[g2].eta2,u72:cases[g72].u72,F:cases[g2].F,combo:cases[g2].combo,nv:false});
    // s_s not entered and the lower bound fails: the seating is unknown, so the
    // station is NOT VERIFIED (blocking) rather than FAIL; its s_s = 0 values are printed
    if(ssDefault && (rec.eta2>1.0001 || rec.u72>1.0001)){
      rec.nv=true;
      rec.msg='NOT VERIFIED: s<sub>s</sub> not entered; at the lower bound s<sub>s</sub> = 0 the station gives F<sub>Ed</sub>/F<sub>Rd</sub> = '+g(rec.eta2,3)+(rec.u72>1.0001? ' and 7.2 = '+g(rec.u72,3) : '');
      out.unsupported.push('Web transverse force at x = '+g(s.x/1000,3)+' m ('+label+'): the stiff bearing length s<sub>s</sub> is not entered; at the lower bound s<sub>s</sub> = 0 the station gives F<sub>Ed</sub>/F<sub>Rd</sub> = '+g(rec.eta2,3)+' (F<sub>Ed</sub> = '+g(rec.F,1)+' kN, F<sub>Rd</sub> = '+g(rec.FRdTot,1)+' kN, type ('+rec.type+'))'+(rec.u72>1.0001? ' and (&eta;<sub>2</sub> + 0.8&eta;<sub>1</sub>)/1.4 = '+g(rec.u72,3) : '')+'. F<sub>Rd</sub> rises with the seating length: enter s<sub>s</sub> (mm along the member, EN 1993-1-5 6.3(1)) in the end row, or tick "bearing stiffener provided"; PASS is blocked until then.');
    }
    out.stations.push(rec);
  });
  const checked=out.stations.filter(s=>!s.stiff && !s.nv);
  out.anyNv=out.stations.some(s=>s.nv);
  out.checked=checked.length>0;
  // station whose derivation is printed: the worst F_Ed/F_Rd of every evaluated
  // station, a NOT VERIFIED one included (its s_s = 0 chain is shown as such)
  const shown=out.stations.filter(s=>!s.stiff);
  if(shown.length){ let w=shown[0]; shown.forEach(s=>{ if(s.eta2>w.eta2) w=s; }); out.show=w; } else out.show=null;
  if(out.checked){
    let w2=checked[0], w72=checked[0];
    checked.forEach(s=>{ if(s.eta2>w2.eta2) w2=s; if(s.u72>w72.u72) w72=s; });
    out.gov2=w2; out.gov72=w72; out.util2=w2.eta2; out.util72=w72.u72;
    checked.filter(s=>s.eta2>1.0001).forEach(s=>out.advisory.push('Bearing stiffener required at x = '+g(s.x/1000,3)+' m ('+s.label+'): F<sub>Ed</sub> = '+g(s.F,1)+' kN exceeds F<sub>Rd</sub> = '+g(s.FRdTot,1)+' kN (EN 1993-1-5 6.2, type ('+s.type+')); tick "bearing stiffener provided" once a stiffener is designed to EN 1993-1-5 9.4, or increase the stiff bearing length s<sub>s</sub>.'));
  }
  return out;
}
/* ---- Shear-reduced cross-section resistances (19 Sep 2026 gap closure, items
   2.10 and 2.13): the yield strength of the shear area is reduced to
   (1 - rho) f_y with rho = (2V_Ed/V_pl(,T),Rd - 1)^2 (cl 6.2.8(3), 6.2.10(3)).
   Returns, in kN.m / kN, the resistances of the section with that reduction:
     MvY   major-axis moment resistance M_v,y,Rd
             I/H Class 1/2 : (W_pl,y - rho A_v^2/(4 t_w)) f_y             (6.2.8(5), Eq 6.30 with A_v for A_w: conservative)
             I/H Class 3   : (W_el,y - rho I_web/(h/2)) f_y                (elastic: web fibre stress limited to (1 - rho) f_y, conservative)
             channel       : (W_pl,y or W_el,y - rho t_w h_w^2/4) f_y     (plastic web modulus; conservative for Class 3)
             RHS/SHS       : (W_pl,y or W_el,y - rho t (h - 2t)^2/2) f_y  (two webs, plastic web modulus; conservative for Class 3)
     MvZ   minor-axis resistance with the web contribution reduced
             I/H           : (W_pl,z - rho A_v t_w/4) f_y (Class 1/2), (W_el,z - rho A_v t_w/6) f_y (Class 3)
             channel       : M_c,z,Rd (1 - rho A_v/A)  (web share of W_z bounded by its area share: conservative)
             RHS/SHS       : (W_pl,z or W_el,z - rho t (h - 2t)(b - t)) f_y (two webs, plastic web modulus about z)
     NV    axial resistance N_V,Rd = (A - rho A_v) f_y/gamma_M0
     aV    the 6.2.9.1 parameter a evaluated on the reduced-yield section:
           ((A - 2 b t_f) - rho A_v)/(A - rho A_v) <= 0.5 (I/H, RHS a_w); a_f,V = (A - 2 h t)/(A - rho A_v) <= 0.5 (RHS)
   Every value is capped at its unreduced counterpart and floored at 0. Pure. */
function shearReducedResistances(sec,cl,fy,Av,rho,gM0){
  gM0=gM0||1.0;
  const A=sec.A*1e2, tw=sec.tw, tf=sec.tf, D=sec.D, B=sec.B, hw=D-2*tf;
  const cls12=cl.cls<=2;
  const Wy=(cls12? sec.Sx : sec.Zx)*1e3, Wz=(cls12? sec.Sy : sec.Zy)*1e3;
  const Mc=Wy*fy/gM0/1e6, Mcz=Wz*fy/gM0/1e6;
  const NV=Math.max(A-rho*Av,0)*fy/gM0/1000;
  let dWy, dWz, form, formZ;
  if(sec.kind==='I'){
    if(cls12){ dWy=rho*Av*Av/(4*tw); dWz=rho*Av*tw/4; form='(W<sub>pl,y</sub> &minus; &rho;A<sub>v</sub>&sup2;/4t<sub>w</sub>)f<sub>y</sub>/&gamma;<sub>M0</sub>'; formZ='(W<sub>pl,z</sub> &minus; &rho;A<sub>v</sub>t<sub>w</sub>/4)f<sub>y</sub>/&gamma;<sub>M0</sub>'; }
    else { const Iweb=tw*Math.pow(hw,3)/12; dWy=rho*Iweb/(D/2); dWz=rho*Av*tw/6; form='(W<sub>el,y</sub> &minus; &rho;I<sub>web</sub>/(h/2))f<sub>y</sub>/&gamma;<sub>M0</sub> [elastic, conservative]'; formZ='(W<sub>el,z</sub> &minus; &rho;A<sub>v</sub>t<sub>w</sub>/6)f<sub>y</sub>/&gamma;<sub>M0</sub>'; }
  } else if(sec.kind==='channel'){
    dWy=rho*tw*hw*hw/4; dWz=rho*Av/A*Wz;
    form='('+(cls12? 'W<sub>pl,y</sub>' : 'W<sub>el,y</sub>')+' &minus; &rho;t<sub>w</sub>h<sub>w</sub>&sup2;/4)f<sub>y</sub>/&gamma;<sub>M0</sub>'+(cls12? '' : ' [plastic web modulus, conservative]');
    formZ='M<sub>c,z,Rd</sub>(1 &minus; &rho;A<sub>v</sub>/A) [web share bounded by its area share, conservative]';
  } else {
    const t=tf, hi=D-2*t;
    dWy=rho*t*hi*hi/2; dWz=rho*t*hi*(B-t);
    form='('+(cls12? 'W<sub>pl,y</sub>' : 'W<sub>el,y</sub>')+' &minus; &rho;t(h &minus; 2t)&sup2;/2)f<sub>y</sub>/&gamma;<sub>M0</sub> [two webs'+(cls12? '' : '; plastic web modulus, conservative')+']';
    formZ='('+(cls12? 'W<sub>pl,z</sub>' : 'W<sub>el,z</sub>')+' &minus; &rho;t(h &minus; 2t)(b &minus; t))f<sub>y</sub>/&gamma;<sub>M0</sub> [two webs]';
  }
  const MvY=Math.min(Math.max((Wy-dWy)*fy/gM0/1e6,0),Mc);
  const MvZ=Math.min(Math.max((Wz-dWz)*fy/gM0/1e6,0),Mcz);
  const Ared=Math.max(A-rho*Av,1e-9);
  const aV=Math.min(Math.max(((A-2*B*tf)-rho*Av)/Ared,0),0.5);              // sec.tf = t for a hollow section
  const afV= sec.isBox? Math.min(Math.max((A-2*D*tf)/Ared,0),0.5) : null;
  return {rho,MvY,MvZ,NV,Mc,Mcz,dWy,dWz,form,formZ,aV,afV,cls12};
}
/* ---- Restraint design forces (19 Sep 2026 gap closure, items 1.7 / 3.17),
   advisory: for every intermediate lateral restraint and every support's
   torsional (fork) restraint, N_f,Ed = M_Ed/h with M_Ed the largest moment at
   the restraint station over the ULS combinations and h the overall depth
   (EN 1993-1-1 5.3.3(3) notation), and the design force 2.5 % N_f,Ed
   (6.3.5.2(5)(b) / SCI practice). Only on the not-fully-restrained path (a
   fully restrained flange has no discrete restraints). Pure. */
function restraintForces(a,sec){
  if((S.restraint||'full')==='full') return null;
  const h=sec.D;
  const pts=[];
  endsList().filter(e=>e.rx).forEach(e=>pts.push({x:e.x*1000,kind:'support',n:e.n,label:'End '+e.n+' (torsional restraint, R<sub>x</sub> held)'}));
  (S.ltbRestraints||[]).forEach((r,i)=>{ const x=(+r.pos)*1000; if(!isFinite(x)||x<-1e-6||x>a.L+1e-6) return; if(r.v===false&&r.phi===false) return;
    pts.push({x,kind:'lateral',n:i+1,label:'lateral restraint '+(i+1)+(r.v===false? ' (twist only)' : '')}); });
  pts.sort((p,q)=>p.x-q.x || (p.kind==='support'? -1 : 1));
  const rows=pts.map(p=>{
    let MEd=0, combo='';
    (a.ulsResults||[]).forEach(res=>{
      const m=Math.max(Math.abs(interpAt(res.fb.xs,res.fb.M,Math.max(p.x-1e-4,0))),Math.abs(interpAt(res.fb.xs,res.fb.M,Math.min(p.x+1e-4,a.L))))/1e6;
      if(m>MEd){ MEd=m; combo=res.combo.label; }
    });
    const NfEd=MEd*1000/h;             // kN: kN.m x 1000 / mm
    return Object.assign({},p,{MEd,combo,NfEd,F:0.025*NfEd});
  });
  const Fmax=rows.reduce((m,r)=>Math.max(m,r.F),0);
  return {h,rows,Fmax,basis:'N<sub>f,Ed</sub> = M<sub>Ed</sub>/h at the restraint station (h = overall depth, EN 1993-1-1 5.3.3(3)); restraint design force 2.5 % N<sub>f,Ed</sub> (6.3.5.2(5)(b), SCI practice) &mdash; advisory, not part of the member verdict; the bracing system must also satisfy the 5.3.3 stiffness/imperfection requirements. A station with M<sub>Ed</sub> = 0 (simply supported end) gets no flange force from this rule; its torsional restraint must still prevent twist.'};
}
/* ===========================================================================
   20 Sep 2026 torsion + N/Mz: combined torsion with direct axial force N_Ed
   and an imposed minor-axis moment M_z,Ed, verified "as Eurocode advises".
   The former block ("Combined torsion with direct axial force or imposed
   minor-axis bending is not implemented as one interaction") is replaced by
   the checks below; torsionCombinedBasis() states the basis in the brief.
   20 Sep 2026 review (reviewer findings on the first implementation): the
   binding policy of (6.1), the closed-section shear flow, the flange shear
   flow at P2, the channel points P1b/P2 and the per-path basis text below.
   ---------------------------------------------------------------------------
   FORMULAE IMPLEMENTED (units: N, mm, N/mm2 unless stated; E = a.E N/mm2,
   G = 81000 N/mm2, gamma_M0 = gamma_M1 = 1.0 UK NA, f_yd = f_y/gamma_M0)
   (1) CROSS-SECTION, EN 1993-1-1 cl 6.2.7(4)-(5) with the yield criterion
       6.2.1(5) Eq (6.1), sigma_z,Ed = 0:
         (sigma_x,Ed/f_yd)^2 + 3 (tau_Ed/f_yd)^2 <= 1
       at every station x of the torsion solution of every ULS combination
       and at the critical points of the section, every contribution taken at
       its worst (|values|, additive - conservative):
         sigma_x,Ed = N_Ed/A + M_y,Ed/W_el,y + M_z,tot/W_el,z + sigma_w,Ed
         tau_Ed     = tau_V + tau_t,Ed + tau_w,Ed
         N_Ed [kN] -> N_Ed*1e3/(A*100)            A = sec.A cm2
         M_z,tot(x) = |M_z,Ed| + |phi(x) M_y,Ed(x)|  imposed S.Mz [kN.m] plus the
                      second-order minor-axis moment of the twist (SCI P385
                      3.1.2 "M_z = phi M_y"; EN 1993-1-1 5.2.1(3): second-order
                      effects included where they increase the action effects)
         sigma_w  = E W_n phi''      W_n = tp.Wn0*100 mm2 (flange tip); channel
                                     tp.Wn2*100 mm2 at the web-flange junction;
                                     I/H: omega = 0 at the junction (web line)
         tau_t    = G t phi'          t = t_f on the flange, t_w on the web
                                     (= T_t t/I_T with T_t = G I_T phi')
         tau_w    = E S_w phi'''/t    S_w = tp.Sw1*1e4 mm4 at the flange centre
                                     (I/H: the junction, where omega = 0);
                                     channel: Sw1 at the flange point where
                                     W_n = 0 (P1b), Sw2 at the junction (flange
                                     side t = t_f, web side t = t_w), Sw3 at the
                                     web mid-depth (P385 Table A.2 points 1/2/3)
         tau_V    = V S/(I_y t)       EN 1993-1-1 6.2.6(4) Eq (6.20), I_y =
                                     sec.Ix*1e4 mm4 (major axis), the first
                                     moment of the area beyond the point:
                                     flange at the junction S_fF = S_f/2 (I/H,
                                     half flange one side of the web) or S_f
                                     (channel, the whole flange), t = t_f;
                                     web at the junction S_f = B t_f (D - t_f)/2,
                                     web mid-depth S_max = S_f + t_w (D/2 - t_f)^2/2,
                                     t = t_w (channel: its own B)
       Points (I/H and channel):
         P1 flange tip:         sigma = N/A + M_y/W_el,y + M_z,tot/W_el,z + E Wn0 phi''
                                tau   = G t_f phi'
         P1b flange, W_n = 0    (channel only) at s_1 = B' Wn0/(Wn0 + Wn2) from
                                the toe, B' = B - t_w/2 (omega is linear along
                                the flange: W_n0 at the toe and W_n2 at the web
                                line have opposite signs, so s_1 = B' - e_0):
                                sigma = N/A + M_y (D/2 - t_f/2)/I_y + M_z,tot |y_toe - s_1|/I_z
                                tau   = V s_1 t_f (D - t_f)/2/(I_y t_f) + G t_f phi' + E Sw1 phi'''/t_f
         P2 junction, flange:   sigma = N/A + M_y (D/2 - t_f/2)/I_y + M_z,tot y_web/I_z + E Wn2 phi''
                                tau   = V S_fF/(I_y t_f) + G t_f phi' + E S_w,J phi'''/t_f
                                (S_w,J = Sw1 for I/H, Sw2 for a channel)
         P3 junction, web:      sigma = N/A + M_y (D/2 - t_f)/I_y + M_z,tot y_web/I_z + E Wn2 phi''
                                tau   = V S_f/(I_y t_w) + G t_w phi' + E Sw2 phi'''/t_w
         P4 web mid-depth:      sigma = N/A + M_z,tot y_web/I_z
                                tau   = V S_max/(I_y t_w) + G t_w phi' + E Sw3 phi'''/t_w
         y_web = distance of the web plane from the minor axis: t_w/2 (I/H);
         channel: the web back c_y = B - I_z/W_el,z (W_el,z tabulated at the
         toe, y_toe = I_z/W_el,z) - the M_z stress on the web (M_z y/I_z),
         omitted by a flange-tip-only summation, is ~44 % of the toe stress for
         a PFC and is kept; the Wn2 and Sw2/Sw3 terms are zero for a doubly
         symmetric I/H.
       Hollow sections, cl 6.2.7(7): warping neglected, tau_t = T_Ed/W_t
         (W_t = tp.Wt*1e3 mm3 or the EN 10210-2 value already used), t = wall;
         the shear flow of V_Ed round the closed mid-line (Eq 6.20, q = 0 at
         the flange mid-width by symmetry): Q_c = ((B - t)/2) t ((D - t)/2) at
         the corner, Q_m = Q_c + t (D - t)^2/8 at the web mid-depth:
         corner:            sigma = N/A + M_y/W_el,y + M_z,tot/W_el,z; tau = tau_t + V Q_c/(I_y t)
         web mid-depth:     sigma = N/A + M_z,tot/W_el,z (the webs are the
                            extreme minor-axis fibres); tau = tau_t + V Q_m/(I_y t)
         flange mid-width:  sigma = N/A + M_y/W_el,y;                  tau = tau_t
         phi(x) of the St Venant solution (a.tors.uls) gives M_z,tot.
       Utilisation "Elastic yield criterion (6.1) with torsion, cl 6.2.7(5)"
       = the largest (6.1) value; reported with the station, point and every
       stress component (tor.elastic). BINDING POLICY (elasticBindingPolicy):
       6.2.7(5) is permissive ("the yield criterion in 6.2.1(5) MAY be
       applied") and 6.2.7(6) permits the plastic moment resistance under
       bending + torsion for Class 1/2 sections with B_Ed from the elastic
       analysis, so (6.1) enters the verdict (c.utils) only where the code
       gives no plastic route: Class 3 sections, and Class 1/2 OPEN sections
       with N_Ed != 0 (no expression of EN 1993-1-1 or P385 combines N_Ed with
       the bimoment plastically). Class 1/2 open sections with N_Ed = 0 are
       verified plastically by the P385 3.1.2 interaction under 6.2.7(6);
       Class 1/2 hollow sections by 6.2.7(7)/(9) (Eq 6.28), 6.2.8(4) and
       6.2.9.1/6.2.10 with V_pl,T,Rd: there (6.1) is computed and printed as
       information (c.info; an advisory when it exceeds 1). The P385 3.1.2
       plastic interaction, V/V_pl,T,Rd (6.2.7(9)), T_Ed/T_Rd and the 6.2.9/
       6.2.10 N-M interactions are kept; 3.1.2 now carries M_z,tot in its
       M_z term.
   (2) MEMBER, EN 1993-6 Annex A Eq (A.1) (informative; stated for I
       sections, applied to channels on the SCI P385 basis), P385 6.2/8.2 form:
         M_y,Ed/M_b,Rd + C_mz M_z,tot/M_z,Rd + k_w k_zw k_alpha M_w,Ed/M_f,Rd <= 1
         k_w = 0.7 - 0.2 M_w,Ed/M_f,Rd, k_zw = 1 - M_z,tot/M_z,Rd,
         k_alpha = 1/(1 - M_y,Ed,max/M_cr); M_z,tot = M_z,Ed + phi M_y per station;
         C_mz = 1.0 whenever an imposed M_z exists (constant diagram, psi = 1,
         Table B.3); otherwise the existing proxy (0.9 / 0.95 / 1.0).
   (3) MEMBER, EN 1993-1-1 cl 6.3.3 Eq 6.61/6.62 (annexB2): M_z,Ed = M_z,Ed
       (imposed) + max over the member of |phi M_y| of the combination's own
       torsion solution; C_mz = 1.0 (Table B.3 upper bound: exact for the
       imposed constant M_z, psi = 1, conservative for the twist-induced
       diagram); everything else unchanged.
   (4) ADVISORY (information only, NOT a utilisation), N_Ed > 0 with torsion:
         N_Ed/(chi_z N_Rk/gM1) + k_zy M_y,Ed/(chi_LT M_y,Rk/gM1)
           + k_zz M_z,tot/(M_z,Rk/gM1) + k_w k_zw k_alpha M_w,Ed/M_f,Rd
       = Eq 6.62 as evaluated (buck.u2, its M_z,Ed already M_z,tot) + the
       largest station value of the (A.1) warping term; k_alpha = 1 when no
       LTB check exists (fully restrained, M_cr -> infinity). "Superposition
       of Eq 6.62 and (A.1) - not a Eurocode expression, information only."
   No Eurocode expression combines N_Ed with warping torsion at member level;
   nothing beyond the clauses named above is invented.
   =========================================================================== */
// 20 Sep 2026: the k_alpha-unbounded state is a FAIL row (LTB governs), printed with the FAIL: prefix
const KALPHA_UNBOUNDED_FAIL='FAIL: M_y,Ed reaches the elastic critical moment M_cr (M_y,Ed >= M_cr): lateral-torsional buckling governs before the torsion interaction can be evaluated - the EN 1993-6 Annex A amplifier k_alpha = 1/(1 - M_y,Ed/M_cr) is unbounded and the utilisation is carried as 99; the member is inadequate as arranged (larger section, shorter unrestrained length or compression-flange restraint).';
const ELASTIC_TORSION_UTIL_NAME="Elastic yield criterion (6.1) with torsion, cl 6.2.7(5)";
const SUPERPOSITION_LABEL="superposition of Eq 6.62 and (A.1) - not a Eurocode expression, information only";
/* 20 Sep 2026 review: is the (6.1) value verdict-binding? (see the comment block, item (1) BINDING POLICY)
   Returns {binding, basis}. */
function elasticBindingPolicy(sec,cl,NEd){
  const cls=cl.cls, hasN=Math.abs(NEd)>1e-9;
  if(cls>=3) return {binding:true, basis:'Class 3 section: no plastic resistance exists, so the elastic verification of EN 1993-1-1 6.2.7(5) with the yield criterion 6.2.1(5) Eq (6.1) is the cross-section check (verdict-binding)'};
  if(!sec.isBox && hasN) return {binding:true, basis:'Class 1/2 open section with N_Ed: EN 1993-1-1 6.2.7(6) admits the plastic resistance for bending + torsion only and no expression of EN 1993-1-1 or SCI P385 combines N_Ed with the bimoment plastically, so the elastic verification 6.2.7(5)/(6.1) is the cross-section check (verdict-binding)'};
  if(sec.isBox) return {binding:false, basis:'Class 1/2 hollow section: the cross-section is verified by the plastic route the code gives - T_Ed/T_Rd (6.2.7(1)/(7), warping neglected), V_Ed/V_pl,T,Rd (6.2.7(9), Eq 6.28), rho from V_pl,T,Rd (6.2.8(4)) and the 6.2.9.1/6.2.10 N-M_y-M_z interaction; the elastic verification of 6.2.7(5) is permissive ("may be applied") and is printed for information'};
  return {binding:false, basis:'Class 1/2 open section without N_Ed: EN 1993-1-1 6.2.7(6) permits the plastic moment resistance under bending + torsion with B_Ed from the elastic analysis (the SCI P385 3.1.2 interaction with M_z,tot); the elastic verification of 6.2.7(5) is permissive ("may be applied") and is printed for information'};
}
/* 20 Sep 2026 review: the basis text of the combined verification, per path (the former single sentence
   named warping torsion and (A.1) for hollow sections and 6.61/6.62 where they are not evaluated).
   o = {box, tension, hasN, hasMz, restrained, buckEvaluated, annexEvaluated, binding} */
function torsionCombinedBasis(o){
  const buckTxt= o.buckEvaluated? ', cl 6.3.3 (6.61/6.62) with the second-order minor-axis moment phi.M_y added to M_z,Ed'
    : o.tension? '; cl 6.3.3 (6.61/6.62) is not evaluated for axial tension (N_t,Rd and the 6.2.9 cross-section interaction apply)'
    : o.restrained? '; cl 6.3.3 (6.61/6.62) does not apply to a fully restrained member without axial compression (the 6.2.9 cross-section interaction with M_z,tot applies)'
    : '';
  if(o.box){
    return 'EN 1993-1-1 6.2.7(7): warping neglected for the hollow section; the verification is 6.2.7(5)/(6.1) at the cross-section with tau_t = T_Ed/W_t, the V_Ed shear flow and N/A + M_y/W_el,y + M_z,tot/W_el,z'+(o.binding? '' : ' (information only: the Class 1/2 plastic route governs)')+', T_Ed/T_Rd (6.2.7(1)), V_Ed/V_pl,T,Rd (6.2.7(9), Eq 6.28), 6.2.8(4)/6.2.9.1/6.2.10 with V_pl,T,Rd'+buckTxt+'; EN 1993-6 (A.1) does not apply to a closed section';
  }
  return 'No expression in EN 1993-1-1 or EN 1993-6 combines N_Ed with warping torsion at member level; the verification is EN 1993-1-1 6.2.7(5)/(6.1) at the cross-section with all stresses'+(o.binding? '' : ' (information only: the Class 1/2 plastic route of 6.2.7(6) / SCI P385 3.1.2 with M_z,tot governs)')+buckTxt+(o.annexEvaluated? ', and EN 1993-6 (A.1) with M_z,Ed = M_z + phi.M_y' : (o.restrained? '; EN 1993-6 (A.1) is not evaluated for a fully restrained member (no lateral-torsional buckling)' : ''));
}
// EN 1993-1-1 6.2.1(5) Eq (6.1) with sigma_z,Ed = 0
function yieldCriterion61(sigmaX,tau,fyd){ return Math.pow(sigmaX/fyd,2)+3*Math.pow(tau/fyd,2); }
// one station's points -> totals and (6.1); pts = [{point, sigmaN, sigmaMy, sigmaMz, sigmaW, tauV, tauT, tauW}]
function elasticPoints61(pts,fyd){
  return pts.map(p=>{ const sigmaX=p.sigmaN+p.sigmaMy+p.sigmaMz+p.sigmaW, tau=p.tauV+p.tauT+p.tauW;
    return Object.assign({},p,{sigmaX,tau,u:yieldCriterion61(sigmaX,tau,fyd)}); });
}
// worst point of a station into the running record (st = {x, combo, ...station data})
function elasticGovern(worst,pts,st){
  let w=worst;
  pts.forEach(p=>{ if(p.u>w.u) w=Object.assign({},st,p,{points:pts}); });
  return w;
}
/* (1) open sections: the P385 Method B solution O = a.torsO (sols per ULS
   combination: sol.xs mm, phi rad, p1 = phi' rad/mm, p2 = phi'' rad/mm2,
   p3 = phi''' rad/mm3; fb.xs/fb.M N.mm/fb.V N coincident). NEd kN (design
   value, |.| taken), MzImp kN.m (imposed |S.Mz|). Returns tor.elastic. */
function torsionElasticOpen(a,sec,fy,gM0,O,NEd,MzImp){
  const E=a.E, G=81000, fyd=fy/gM0;
  const A=sec.A*100, Iy=sec.Ix*1e4, Iz=sec.Iy*1e4, Wely=sec.Zx*1e3, Welz=sec.Zy*1e3;   // mm2, mm4, mm3
  const D=sec.D, B=sec.B, tf=sec.tf, tw=sec.tw, chan=sec.kind==='channel', tp=sec.tp||{};
  const Wn0=(tp.Wn0||0)*100, Wn2=chan? (tp.Wn2||0)*100 : 0;                       // mm2 (cm2 -> mm2)
  const Sw1=(tp.Sw1||0)*1e4, Sw2=chan? (tp.Sw2||0)*1e4 : 0, Sw3=chan? (tp.Sw3||0)*1e4 : 0;   // mm4 (cm4 -> mm4)
  const SwJ=chan? Sw2 : Sw1;                                                       // junction, flange side (20 Sep 2026 review: Sw2 for a channel, not max(Sw1,Sw2,Sw3))
  const Sf=B*tf*(D-tf)/2, Smax=Sf+tw*Math.pow(D/2-tf,2)/2;                        // mm3, first moments about y-y (web values, Eq 6.20)
  const SfF=chan? Sf : Sf/2;                                                       // mm3, flange at the junction: one side of the web (20 Sep 2026 review)
  const yWeb=chan? Math.max(B-Iz/Welz,0) : tw/2;                                   // mm, web plane from the minor axis
  // channel point P1b: the flange point where W_n = 0 (P385 point 1, S_w1), s_1 from the toe (omega linear along the flange)
  const Bf=B-tw/2, s1=(chan && Wn0+Wn2>0)? Bf*Wn0/(Wn0+Wn2) : 0, yToe=Iz/Welz, y1b=Math.abs(yToe-s1);
  const sN=Math.abs(NEd)*1e3/A;
  let worst={u:-1}, MzTwistMax=0, MzTotMax=0, nSt=0;
  O.sols.forEach(se=>{ const g=se.sol, fb=se.fb;
    g.xs.forEach((x,i)=>{
      const MyN=Math.abs(interpAt(fb.xs,fb.M,x)), VN=Math.abs(interpAt(fb.xs,fb.V,x));   // N.mm, N
      const phi=g.phi[i], p1=Math.abs(g.p1[i]), p2=Math.abs(g.p2[i]), p3=Math.abs(g.p3[i]);
      const MzTw=Math.abs(phi*MyN), MzTot=MzImp*1e6+MzTw;                          // N.mm
      const sMyTip=MyN/Wely, sMyFm=MyN*(D/2-tf/2)/Iy, sMyJ=MyN*(D/2-tf)/Iy;
      const sMzTip=MzTot/Welz, sMzWeb=MzTot*yWeb/Iz, sMz1b=MzTot*y1b/Iz;
      const sW0=E*Wn0*p2, sW2=E*Wn2*p2;
      const tTf=G*tf*p1, tTw=G*tw*p1;
      const tWJ=E*SwJ*p3/tf, tW1b=chan? E*Sw1*p3/tf : 0, tW3=chan? E*Sw2*p3/tw : 0, tW4=chan? E*Sw3*p3/tw : 0;
      const tV2=VN*SfF/(Iy*tf), tV1b=chan? VN*s1*tf*(D-tf)/2/(Iy*tf) : 0, tV3=VN*Sf/(Iy*tw), tV4=VN*Smax/(Iy*tw);
      const list=[
        {point:'P1 flange tip',                 sigmaN:sN, sigmaMy:sMyTip, sigmaMz:sMzTip, sigmaW:sW0, tauV:0,    tauT:tTf, tauW:0}];
      if(chan) list.push({point:'P1b flange at W_n = 0',   sigmaN:sN, sigmaMy:sMyFm,  sigmaMz:sMz1b,  sigmaW:0,   tauV:tV1b, tauT:tTf, tauW:tW1b});
      list.push(
        {point:'P2 web-flange junction, flange',sigmaN:sN, sigmaMy:sMyFm,  sigmaMz:sMzWeb, sigmaW:sW2, tauV:tV2,  tauT:tTf, tauW:tWJ},
        {point:'P3 web-flange junction, web',   sigmaN:sN, sigmaMy:sMyJ,   sigmaMz:sMzWeb, sigmaW:sW2, tauV:tV3,  tauT:tTw, tauW:tW3},
        {point:'P4 web mid-depth',              sigmaN:sN, sigmaMy:0,      sigmaMz:sMzWeb, sigmaW:0,   tauV:tV4,  tauT:tTw, tauW:tW4});
      const pts=elasticPoints61(list,fyd);
      worst=elasticGovern(worst,pts,{x,combo:se.combo.label,My:MyN/1e6,V:VN/1e3,phi,p1:g.p1[i],p2:g.p2[i],p3:g.p3[i],MzImp,MzTwist:MzTw/1e6,MzTot:MzTot/1e6});
      MzTwistMax=Math.max(MzTwistMax,MzTw/1e6); MzTotMax=Math.max(MzTotMax,MzTot/1e6); nSt++;
    });
  });
  if(worst.u<0) return null;
  return Object.assign(worst,{box:false,fy:fyd,NEd,MzTwistMax,MzTotMax,nStations:nSt,
    geom:{A,Iy,Iz,Wely,Welz,Wn0,Wn2,Sw1,Sw2,Sw3,SwJ,Sf,SfF,Smax,yWeb,yToe,s1,y1b,E,G,tf,tw,D,B,chan},
    basis:'EN 1993-1-1 6.2.7(5): elastic verification of the cross-section with the yield criterion 6.2.1(5) Eq (6.1), sigma_x = N/A + M_y/W_el,y + M_z,tot/W_el,z + sigma_w and tau = tau_V + tau_t + tau_w at the flange tip, '+(chan? 'the flange point where W_n = 0 (S_w1), ' : '')+'the web-flange junction (flange and web side) and the web mid-depth, every station of every ULS combination, every contribution at its worst; M_z,tot = M_z,Ed + phi.M_y'});
}
/* (1) hollow sections, cl 6.2.7(7): warping neglected; the St Venant solution
   a.tors.uls per ULS combination (r.xs/r.T N.mm torque diagram, r.nodes/r.phi
   twist); Wt mm3. The V_Ed shear flow round the closed mid-line (Eq 6.20):
   Q_c at the corner, Q_m at the web mid-depth (20 Sep 2026 review: the corner
   carried tau_V = 0 and the web the mean V/(2(D - 2t)t)). Returns tor.elastic. */
function torsionElasticBox(a,sec,fy,gM0,Wt,NEd,MzImp){
  const fyd=fy/gM0, A=sec.A*100, Iy=sec.Ix*1e4, Wely=sec.Zx*1e3, Welz=sec.Zy*1e3, t=sec.tf, D=sec.D, B=sec.B;
  const Qc=((B-t)/2)*t*((D-t)/2), Qm=Qc+t*Math.pow(D-t,2)/8;                       // mm3, mid-line first moments (cut at the flange mid-width, q = 0)
  const sN=Math.abs(NEd)*1e3/A;
  const torsUls=(a.tors&&a.tors.uls)||[];
  let worst={u:-1}, MzTwistMax=0, MzTotMax=0, nSt=0;
  a.ulsResults.forEach(res=>{
    const tr=torsUls.find(u=>u.combo===res.combo); if(!tr) return;
    const xs=[...new Set([...res.fb.xs,...tr.r.xs].map(x=>+x.toFixed(4)))].sort((p,q)=>p-q);
    xs.forEach(x=>{
      const MyN=Math.abs(interpAt(res.fb.xs,res.fb.M,x)), VN=Math.abs(interpAt(res.fb.xs,res.fb.V,x));
      const T=Math.abs(interpAt(tr.r.xs,tr.r.T,x)), phi=interpAt(tr.r.nodes,tr.r.phi,x);
      const MzTw=Math.abs(phi*MyN), MzTot=MzImp*1e6+MzTw;
      const tauT=T/Wt, tauVc=VN*Qc/(Iy*t), tauVm=VN*Qm/(Iy*t), sMy=MyN/Wely, sMz=MzTot/Welz;
      const pts=elasticPoints61([
        {point:'corner',           sigmaN:sN, sigmaMy:sMy, sigmaMz:sMz, sigmaW:0, tauV:tauVc, tauT, tauW:0},
        {point:'web mid-depth',    sigmaN:sN, sigmaMy:0,   sigmaMz:sMz, sigmaW:0, tauV:tauVm, tauT, tauW:0},
        {point:'flange mid-width', sigmaN:sN, sigmaMy:sMy, sigmaMz:0,   sigmaW:0, tauV:0,     tauT, tauW:0}],fyd);
      worst=elasticGovern(worst,pts,{x,combo:res.combo.label,My:MyN/1e6,V:VN/1e3,T:T/1e6,phi,MzImp,MzTwist:MzTw/1e6,MzTot:MzTot/1e6});
      MzTwistMax=Math.max(MzTwistMax,MzTw/1e6); MzTotMax=Math.max(MzTotMax,MzTot/1e6); nSt++;
    });
  });
  if(worst.u<0) return null;
  return Object.assign(worst,{box:true,fy:fyd,NEd,MzTwistMax,MzTotMax,nStations:nSt,geom:{A,Iy,Wely,Welz,Wt,t,D,B,Qc,Qm},
    basis:'EN 1993-1-1 6.2.7(5) with 6.2.7(7) (warping neglected, tau_t = T_Ed/W_t): yield criterion 6.2.1(5) Eq (6.1) at the corner (tau_t + V Q_c/(I_y t)), the web mid-depth (tau_t + V Q_m/(I_y t)) and the flange mid-width, the V_Ed shear flow round the closed mid-line per Eq (6.20), every station of every ULS combination, every contribution at its worst; M_z,tot = M_z,Ed + phi.M_y'});
}
/* (3) max over the member of |phi(x) M_y,Ed(x)| (kN.m) for ONE ULS combination:
   the P385 solution for open sections, the St Venant twist for hollow ones. */
function torsionMzTwistMax(a,combo){
  let m=0;
  if(a.torsO&&a.torsO.ok&&a.torsO.sols){
    const se=a.torsO.sols.find(s=>s.combo===combo);
    if(se) se.sol.xs.forEach((x,i)=>{ m=Math.max(m,Math.abs(se.sol.phi[i]*interpAt(se.fb.xs,se.fb.M,x))/1e6); });
  } else if(a.tors&&a.tors.on&&a.tors.uls){
    const tr=a.tors.uls.find(u=>u.combo===combo), res=(a.ulsResults||[]).find(r=>r.combo===combo);
    if(tr&&res) tr.r.nodes.forEach((x,i)=>{ m=Math.max(m,Math.abs(tr.r.phi[i]*interpAt(res.fb.xs,res.fb.M,x))/1e6); });
  }
  return m;
}
/* (2) EN 1993-6 Annex A Eq (A.1) over the P385 station grids of tor (rows
   carry My, Mw, Mz = phi.My, MzImp, MzTot kN.m); MbA = chi_LT W_y f_y/gM1
   (no f-factor, P385 basis), McrA kN.m, CmzProxy the caller's diagram proxy.
   Shared by the standard (closed-form Mcr) and eigen paths. */
function annexAEval(tor,MbA,McrA,CmzProxy){
  const MzImp=tor.MzImp||0;
  const Cmz= MzImp>1e-9? 1.0 : CmzProxy;
  const CmzBasis= MzImp>1e-9? 'C_mz = 1.0: the imposed M_z is a constant diagram (psi = 1, Table B.3); the twist-induced part carries no proxy'
                            : 'C_mz = '+CmzProxy.toFixed(2)+' for the twist-induced M_z = phi.M_y diagram (no M_z imposed; the caller\'s load-case proxy or override)';
  let MyMax=0; tor.grids.forEach(g2=>g2.rows.forEach(r2=>{ MyMax=Math.max(MyMax,r2.My); }));
  // resistances CLASS-CONSISTENT (elastic for Class 3)
  const MzR=tor.cls12? tor.Mplz : tor.Melz;
  const MfR=tor.cls12? tor.Mplf : tor.Melf;
  if(MyMax>=McrA*0.999) return {u:99,kAlpha:Infinity,Cmz,CmzBasis,MbA,McrA,MzR,MfR,MzImp,MzTwist:tor.MzMax,MzTot:tor.MzTot,unbounded:true};
  const kAlpha=1/(1-MyMax/McrA);
  let worst={u:-1};
  tor.grids.forEach(g2=>g2.rows.forEach(r2=>{
    const MzT=r2.MzTot;                                     // M_z,tot = M_z,Ed + phi.M_y at the station
    const kw=Math.max(0.7-0.2*r2.Mw/MfR,0);
    const kzw=Math.max(1-MzT/MzR,0);
    const u=r2.My/MbA + Cmz*MzT/MzR + kw*kzw*kAlpha*r2.Mw/MfR;
    if(u>worst.u) worst={u,x:r2.x,My:r2.My,Mz:MzT,MzImp,MzTwist:r2.Mz,MzTot:MzT,Mw:r2.Mw,kw,kzw,combo:g2.combo.label};
  }));
  return Object.assign(worst,{kAlpha,Cmz,CmzBasis,MbA,McrA,MzR,MfR,unbounded:false});
}
/* (4) advisory superposition, N_Ed > 0 with an open-section torsion solution;
   annex = null when no LTB check exists (k_alpha = 1). Returns null or
   {u, u62, uw, kAlpha, at, text}; the caller pushes text as an advisory. */
function torsionSuperposition(tor,annex,buck){
  if(!(tor&&tor.p385&&tor.grids&&buck&&buck.Fc>1e-9)) return null;
  if(annex&&annex.unbounded) return null;
  const MzR=tor.cls12? tor.Mplz : tor.Melz, MfR=tor.cls12? tor.Mplf : tor.Melf;
  const kAlpha=(annex&&isFinite(annex.kAlpha))? annex.kAlpha : 1;
  let uw=0, at=null;
  tor.grids.forEach(g2=>g2.rows.forEach(r2=>{
    const kw=Math.max(0.7-0.2*r2.Mw/MfR,0), kzw=Math.max(1-r2.MzTot/MzR,0);
    const t=kw*kzw*kAlpha*r2.Mw/MfR;
    if(t>uw){ uw=t; at={x:r2.x,kw,kzw,Mw:r2.Mw,MzTot:r2.MzTot,combo:g2.combo.label}; }
  }));
  const u=buck.u2+uw;
  const text='ADVISORY - '+SUPERPOSITION_LABEL+': N<sub>Ed</sub>/(&chi;<sub>z</sub>N<sub>Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>zy</sub>M<sub>y,Ed</sub>/(&chi;<sub>LT</sub>M<sub>y,Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>zz</sub>M<sub>z,tot</sub>/(M<sub>z,Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>w</sub>k<sub>zw</sub>k<sub>&alpha;</sub>M<sub>w,Ed</sub>/M<sub>f,Rd</sub> = '
    +buck.u2.toFixed(3)+' (Eq 6.62 with M<sub>z,tot</sub> = '+buck.MzEd.toFixed(2)+' kN&middot;m) + '+uw.toFixed(3)+' (warping term of (A.1)'+(at? ' at x = '+(at.x/1000).toFixed(2)+' m, k<sub>w</sub> = '+at.kw.toFixed(3)+', k<sub>zw</sub> = '+at.kzw.toFixed(3)+', k<sub>&alpha;</sub> = '+kAlpha.toFixed(3)+(annex? '' : ' (no LTB check: M<sub>cr</sub> unbounded)') : '')+') = '+u.toFixed(3)+'. Not a utilisation; the verdict rests on (6.1), 6.61/6.62, (A.1), V/V<sub>pl,T,Rd</sub> and P385 3.1.2.';
  return {u,u62:buck.u2,uw,kAlpha,at,text,label:SUPERPOSITION_LABEL};
}
function checksEC3Restrained(a){
  // SCI worked-example procedure: fully laterally restrained beam to BS EN 1993-1-1 (UK NA).
  // Sequence: classification -> shear (6.2.6) -> shear buckling screen (6.2.6(6)) ->
  // moment (6.2.5, with the 6.2.8 shear check made AT the point of maximum moment) ->
  // vertical deflection (NA 2.23).
  const sec=a.sec, fy=a.py;
  const unsupported=[];
  const advisory=[];
  const eps=epsEC3(fy);
  const gM0=1.0, eta=1.0; // UK NA gammaM0; eta=1.0 taken conservatively (EN 1993-1-5)
  const F=S.axial||0;
  // classification must see the coexistent axial compression: with N_Ed the web
  // is an "internal part in bending and compression" (Table 5.2) and its limits
  // tighten below 72e/83e/124e - a web that is Class 1 in pure bending can be
  // Class 3 or 4 under combined actions, which changes W_y and can invalidate
  // the plastic M_N,Rd expressions.
  const cl=classifyEC3(sec,eps,{NEd:Math.max(F,0)*1000, fy,minorBending:Math.abs(S.Mz||0)>1e-9});
  // I/H under M_z (item 1.9(b)): the web is unstressed by M_z and keeps its y-y
  // (+N) limits; the flange outstands are classified under the combined stress
  // (the outstand compressed by M_z is wholly in compression, alpha = 1, so the
  // 9e/10e/14e bound governs). Printed with the classification.
  cl.mzStress = (Math.abs(S.Mz||0)>1e-9 && sec.kind==='I')? mzFlangeStress(sec,F,a.Mmax,S.Mz) : null;
  if(Math.abs(S.Mz||0)>1e-9 && sec.kind!=='I') advisory.push('For biaxial bending of this section family each wall parallel to the web is conservatively classified using the uniform-compression limits in Table 5.2; no favourable biaxial stress distribution is assumed.');
  // Class-4 web in uniform compression (item 1.9(c) / 3.9(c)): effective area
  // per EN 1993-1-5 4.4 (aeffWebCompression) instead of blocking; used in
  // N_c,Rd, N_b,Rd and the Table 6.7 Class-4 column of Eq 6.61/6.62. The
  // flanges must stay Class <= 3 in uniform compression and the section must
  // be doubly symmetric (e_N = 0); a channel web reduction would shift the
  // centroid (N e_N), which is not modelled.
  const aeff=aeffWebCompression(sec,eps);
  aeff.active = F>0 && aeff.applies;
  const flangeCompLimit = sec.isBox? 42*eps : 14*eps;
  if(F>0 && sec.bT>flangeCompLimit) unsupported.push('The flange is Class 4 in uniform compression (c/t = '+g(sec.bT,2)+' > '+g(flangeCompLimit,2)+'): an effective flange width per EN 1993-1-5 4.4 is not implemented; PASS is blocked.');
  if(aeff.active && sec.kind==='channel') unsupported.push('PFC web Class 4 in uniform compression (d/t = '+g(sec.dt,2)+' > 42&epsilon; = '+g(42*eps,2)+'): the effective web shifts the centroid of a channel (N<sub>Ed</sub>e<sub>N</sub> minor-axis moment, EN 1993-1-1 6.2.9.3), which is not implemented; PASS is blocked.');
  const clsName=["","Class 1","Class 2","Class 3","Class 4"][cl.cls];
  if(cl.cls>=4) unsupported.push("EC3 Class 4 (slender) section"+(cl.webCase==='bending+compression'?" (web classified for combined bending + compression: the effective section under the actual stress gradient, with its e_N shift, is not implemented)":"")+": effective-section properties per EN 1993-1-5 are required; not covered by the restrained-beam procedure.");
  // web transverse forces at every point load and every support reaction
  // (EN 1993-1-5 clause 6 + 7.2; webTransverseCheck above); its utilisations
  // enter the verdict below, a declared stiffener prints an advisory
  const web=webTransverseCheck(a,sec,fy,eps,cl);
  web.unsupported.forEach(m=>unsupported.push(m));
  web.advisory.forEach(m=>advisory.push(m));
  advisory.push("Web transverse forces (EN 1993-1-5 clause 6) are checked at every point load and support as patch loads on the loaded flange: distributed loads, loads hung from the bottom flange (hangers), the total-load check of closely spaced loads (6.3(3)) and flange-induced buckling (section 8) are not evaluated; a declared bearing stiffener must be designed to 9.4.");
  const MzEd=Math.abs(S.Mz||0);   // applied minor-axis design moment (kN.m), single-value input
  // ---- axial + biaxial bending cross-section resistance, cl 6.2.9.1
  // Forms (My/MN,y)^alpha + (Mz/MN,z)^beta <= 1. (validated against the
  // independent commercial-software SHS beam-column example for the uniaxial+axial case) ----
  const AgAx=sec.A*1e2, AnetAx=(S.anet!=null? S.anet*1e2 : AgAx);
  const fuAx=fuFromGrade(S.grade);
  const NplRd=AgAx*fy/gM0/1000;                              // kN
  const NcRd= aeff.active? aeff.Aeff*fy/gM0/1000 : NplRd;    // kN, cl 6.2.4(2): A_eff for a Class-4 web in uniform compression
  const NuRd=0.9*AnetAx*fuAx/1.10/1000;                      // kN, gammaM2 = 1.10 (UK NA Table NA.1: resistance of cross-sections in tension to fracture; 1.25 is the EN recommended value)
  const NtRd=Math.min(NplRd,NuRd);
  let ax=null;
  if(Math.abs(F)>1e-9 || MzEd>1e-9){
    const nAx=Math.abs(F)/NplRd;
    const biax=MzEd>1e-9;
    const My=Math.abs(a.Mmax);
    // minor-axis resistances reported with the axial/biaxial block (so the
    // brief prints engine values, not its own): V_pl,z,Rd from the cl 6.2.6(3)
    // shear area for load parallel to the flanges - (f) rolled I, H and
    // channel sections A - sum(hw tw), (h) hollow sections A b/(b + h) - and
    // M_c,z,Rd = W_z fy/gammaM0 with the class-consistent modulus (cl 6.2.5).
    // V_z,Ed is zero in this single-plane model, so V_pl,z,Rd is informational.
    const hwAxz=sec.D-2*sec.tf;
    const Avz= sec.isBox? AgAx*sec.B/(sec.D+sec.B) : Math.max(AgAx-hwAxz*sec.tw,0);
    const VplZ=Avz*fy/(Math.sqrt(3)*gM0)/1000;                 // kN
    const Mcz=(cl.cls<=2? sec.Sy : sec.Zy)*1e3*fy/gM0/1e6;      // kN.m
    if(sec.kind==='channel' && cl.cls<4){
      const McChan=(cl.cls<=2? sec.Sx:sec.Zx)*1e3*fy/gM0/1e6;
      const McChanZ=Mcz;
      ax={cls3:false,chan:true,n:nAx,NplRd,NuRd,NtRd,MN:McChan,MNz:McChanZ,alpha:1,beta:1,biax,Mz:MzEd,Avz,VplZ,Mcz,
        mnLbl:'channel: linear interaction (cl 6.2.1(7); &alpha; = &beta; = 1, conservative &mdash; verified commercial-software basis)',
        mUtil:nAx + My/Math.max(McChan,1e-9) + (biax? MzEd/Math.max(McChanZ,1e-9):0),
        nUtil: F>0? Math.abs(F)/NcRd : (F<0? Math.abs(F)/NtRd : 0), tension:F<0, NcRd, aeff};
    } else if(cl.cls>=4){
      /* Class 4 already blocked above */
    } else {
      const Mpl=sec.Sx*1e3*fy/gM0/1e6, Mel=sec.Zx*1e3*fy/gM0/1e6;
      const Mplz=sec.Sy*1e3*fy/gM0/1e6, Melz=sec.Zy*1e3*fy/gM0/1e6;
      let MN,MNz,alpha,beta,mnLbl,waiver=false;
      if(cl.cls<=2){
        const aw=Math.min(Math.max((AgAx-2*sec.B*sec.tf)/AgAx,0),0.5);
        MN=Math.max(0,Math.min(Mpl*(1-nAx)/(1-0.5*aw),Mpl));
        if(sec.kind==='I'){
          const hwAx=sec.D-2*sec.tf;
          if(Math.abs(F)*1000<=0.25*NplRd*1000 && Math.abs(F)*1000<=0.5*hwAx*sec.tw*fy/gM0){ MN=Mpl; waiver=true; }
          alpha=2; beta=Math.max(5*nAx,1);
          MNz = nAx<=aw ? Mplz : Mplz*(1-Math.pow((nAx-aw)/(1-aw),2));  // cl 6.2.9.1(5), minor axis
          mnLbl='a = '+g(aw,3)+'; &alpha; = 2, &beta; = '+g(beta,2)+(waiver? '; small axial (cl 6.2.9.1(4)): no major reduction':'');
        } else { // RHS / SHS
          alpha=nAx<=0.8?1.66/(1-1.13*nAx*nAx):6; beta=alpha;
          const af=Math.min(Math.max((AgAx-2*sec.D*sec.tf)/AgAx,0),0.5);
          MNz=Math.min(Mplz*(1-nAx)/(1-0.5*af),Mplz);
          mnLbl='a<sub>w</sub> = '+g(aw,3)+'; &alpha; = &beta; = 1.66/(1&minus;1.13n&sup2;) = '+g(alpha,2);
        }
        ax={cls3:false,n:nAx,NplRd,NuRd,NtRd,MN,MNz,alpha,beta,mnLbl,biax,Mz:MzEd,Avz,VplZ,Mcz,
          mUtil:Math.pow(My/Math.max(MN,1e-9),alpha) + (biax? Math.pow(MzEd/Math.max(MNz,1e-9),beta):0),
          nUtil: F>0? Math.abs(F)/NcRd : (F<0? Math.abs(F)/NtRd : 0), tension:F<0, NcRd, aeff};
      } else {
        ax={cls3:true,n:nAx,NplRd,NuRd,NtRd,MN:Mel,MNz:Melz,alpha:1,beta:1,mnLbl:'Class 3: elastic, cl 6.2.9.2',biax,Mz:MzEd,Avz,VplZ,Mcz,
          mUtil:Math.abs(F)/NplRd + My/Math.max(Mel,1e-9) + (biax? MzEd/Math.max(Melz,1e-9):0),
          nUtil: F>0? Math.abs(F)/NcRd : (F<0? Math.abs(F)/NtRd : 0), tension:F<0, NcRd, aeff};
      }
    }
  }
  // geometry
  const hw=sec.D-2*sec.tf;                                   // clear web depth h - 2tf
  const cOut = sec.kind==='I'? (sec.B-sec.tw-2*sec.r)/2 : null; // outstand flange width
  // shear area, cl 6.2.6(3)
  const A=sec.A*100; let Av,AvRaw=null,avFloor=null;
  if(sec.isBox){ Av=A*sec.D/(sec.D+sec.B); }
  else if(sec.kind==='channel'){ Av=A-2*sec.B*sec.tf+(sec.tw+sec.r)*sec.tf; }
  else { AvRaw=A-2*sec.B*sec.tf+(sec.tw+2*sec.r)*sec.tf; avFloor=eta*hw*sec.tw; Av=Math.max(AvRaw,avFloor); }
  const VcRd=Av*fy/(Math.sqrt(3)*gM0)/1000;                   // kN
  const Fv=Math.abs(a.Vmax);
  const shearUtil=VcRd>0? Fv/VcRd : 0;
  const cshear=(fy/Math.sqrt(3))/gM0;
  const shearTorsionUtil=(VplTRd,V=Fv)=> VplTRd>1e-9 ? V/VplTRd : (V>1e-9 ? 99 : 0);
  const sweepXs=(...sets)=>[...new Set(sets.flat().map(x=>+x.toFixed(4)))].sort((p,q)=>p-q);
  const boxShearTorsionSweep=(Wt)=>{
    let worst={u:-1,x:0,V:0,T:0,tau:0,VplTRd:VcRd,combo:''};
    const torsUls=(a.tors&&a.tors.uls)||[];
    a.ulsResults.forEach(res=>{
      const tr=torsUls.find(u=>u.combo===res.combo);
      if(!tr) return;
      sweepXs(res.fb.xs,tr.r.xs).forEach(x=>{
        const Vx=Math.abs(interpAt(res.fb.xs,res.fb.V,x))/1000;
        const T=Math.abs(interpAt(tr.r.xs,tr.r.T,x))/1e6;
        const tau=T*1e6/Wt;
        const VplTRd=Math.max(0,1-tau/cshear)*VcRd;
        const u=shearTorsionUtil(VplTRd,Vx);
        if(u>worst.u) worst={u,x,V:Vx,T,tau,VplTRd,combo:res.combo.label,zeroCapacity:VplTRd<=1e-9&&Vx>1e-9};
      });
    });
    if(worst.u<0){
      const T=(a.tors&&a.tors.Tmax)||0, tau=T*1e6/Wt;
      const VplTRd=Math.max(0,1-tau/cshear)*VcRd;
      worst={u:shearTorsionUtil(VplTRd,Fv),x:(a.tors&&a.tors.Tpos? a.tors.Tpos:0)*1000,V:Fv,T,tau,VplTRd,combo:(a.tors&&a.tors.governT)||'',zeroCapacity:VplTRd<=1e-9&&Fv>1e-9};
    }
    return worst;
  };
  // shear buckling screen, cl 6.2.6(6): hw/tw <= 72 eps/eta (unstiffened web)
  const sbRatio = sec.isBox? sec.dt : hw/sec.tw;
  const sbLimit = 72*eps/eta;
  const sbOk = sbRatio<=sbLimit+1e-9;
  if(!sbOk) unsupported.push("h_w/t_w exceeds 72*eps/eta: web shear-buckling resistance per EN 1993-1-5 must be checked and is not implemented in this calculator.");
  // torsion (cl 6.2.7) - active when loads are eccentric to the shear centre
  let tor=null;
  if(a.tors && a.tors.on){
    if(sec.isBox){
      const cb=ctBoxEN10210(sec);
      const Wt=(sec.tp&&sec.tp.Wt)? sec.tp.Wt*1e3 : cb.Ct; // mm3, tabulated P385 value preferred
      const WtSrc=(sec.tp&&sec.tp.Wt)? 'SCI P385 Table A.7/A.8' : 'EN 10210-2 formula';
      const ItShow=(sec.tp&&sec.tp.IT)? sec.tp.IT*1e4 : cb.It;
      const TRd=fy*Wt/(Math.sqrt(3)*gM0)/1e6;                  // kN.m, cl 6.2.7(7); Tw,Ed neglected for hollow sections
      const TEd=a.tors.Tmax;
      const tauMax=TEd*1e6/Wt;                                  // N/mm2, shear stress due to peak torsion
      const vt=boxShearTorsionSweep(Wt);
      // 20 Sep 2026 torsion + N/Mz: elastic yield criterion (6.1) with N_Ed, M_y, M_z,tot and
      // tau_t = T_Ed/W_t (cl 6.2.7(5)/(7)), every station of every ULS combination (helper above)
      const elastic=torsionElasticBox(a,sec,fy,gM0,Wt,F,MzEd);
      tor={box:true,Wt,WtSrc,ItShow,TRd,TEd,
        tau:vt.tau, tauMax, VplTRd:vt.VplTRd, vt,
        torUtil:TRd>0? TEd/TRd:0, vtUtil:vt.u, vtZeroCapacity:!!vt.zeroCapacity,
        elastic, MzImp:MzEd, MzTot:elastic? elastic.MzTotMax : MzEd, MzTwistMax:elastic? elastic.MzTwistMax : 0,
        GIt:a.tors.GIt,TmaxSLS:a.tors.TmaxSLS,phiMax:a.tors.phiMax,phiDeg:a.tors.phiMax*180/Math.PI,phiPos:a.tors.phiPos,governT:a.tors.governT,governTw:a.tors.governTw};
      if(sec.boxType==='CF') unsupported.push("Torsional constants are computed with hot-finished (EN 10210-2) corner geometry; cold-formed (EN 10219-2) corners differ slightly - verify W_t for a cold-formed section.");
    } else if(a.torsO && a.torsO.ok){
      // ---- SCI P385 Method B: elastic warping analysis + design effects ----
      const O=a.torsO, hh=sec.D-sec.tf;                     // (h - tf) flange lever
      const EIw=a.E*O.Iw;
      const chan=sec.kind==='channel';
      const Mply=sec.Sx*1e3*fy/1e6, Mplz=sec.Sy*1e3*fy/1e6; // kNm
      const Mely=sec.Zx*1e3*fy/1e6, Melz=sec.Zy*1e3*fy/1e6;
      // 20 Sep 2026 accuracy: M_f,Rd = the resistance of ONE FLANGE bending about its own axis in the flange
      // plane (SCI P385 3.1.2 / 6.2, W_pl,f = t_f b^2/4, W_el,f = t_f b^2/6) for I/H AND channels; formerly
      // M_pl,z/2 for I/H, which carried half the web's minor-axis share (UB 457x191x82: 41.80 vs 40.255 kN.m,
      // reference software prints 40.255) - affects P385 3.1.2, the (A.1) warping term and k_w by up to 4 %
      const Mplf=sec.B*sec.B*sec.tf/4*fy/1e6;               // one flange (plastic), kN.m
      const Melf=sec.B*sec.B*sec.tf/6*fy/1e6;               // one flange (elastic), kN.m
      const cls12=cl.cls<=2;
      // evaluate effects on each ULS combo's own coincident (My, phi, Mw) fields
      const SwChan=(chan&&sec.tp)? Math.max(sec.tp.Sw2||0,sec.tp.Sw3||0)*1e4 : 0;
      const MzImp=MzEd;                                       // 20 Sep 2026 torsion + N/Mz: imposed constant M_z,Ed (kN.m, S.Mz)
      let cross={u:-1}, grids=[], tauT=0, tauW=0, TtEnds=[0,0], TEnds=[0,0], BMaxAbs=0, BMaxPos=0, MwMaxAbs=0, MzMax=0, MzTotMax=0, phiUmax=0, MyAtCross=0;
      let vt={u:-1,x:0,V:0,T:0,tauT:0,tauW:0,VplTRd:VcRd,combo:'',zeroCapacity:false};
      O.sols.forEach(se=>{
        const g=se.sol, fb=se.fb;
        const rows=g.xs.map((x,i)=>{
          const My=Math.abs(interpAt(fb.xs,fb.M,x))/1e6;         // kNm
          const Vx=Math.abs(interpAt(fb.xs,fb.V,x))/1000;        // kN, coincident shear
          const phi=g.phi[i];
          const Mw=Math.abs(EIw*g.p2[i]/hh)/1e6;                 // kNm, per flange
          const Mz=Math.abs(phi*interpAt(fb.xs,fb.M,x)/1e6);     // kNm (phi*My), second-order minor-axis moment of the twist
          const MzTot=MzImp+Mz;                                  // 20 Sep 2026 torsion + N/Mz: M_z,tot = M_z,Ed + phi.M_y (P385 3.1.2 M_z term; = Mz when no M_z is imposed)
          const Tt=Math.abs(O.GIt*g.p1[i])/1e6;                  // kNm, coincident St Venant torque
          const tauTi=Tt*1e6*sec.tw/O.IT;                        // N/mm2
          const tauWi=SwChan? Math.abs(a.E*SwChan*g.p3[i]/sec.tw) : 0;
          const VplTRdi=chan
            ? Math.max(0,(Math.sqrt(Math.max(0,1-tauTi/(1.25*cshear)))-tauWi/cshear))*VcRd
            : Math.sqrt(Math.max(0,1-tauTi/(1.25*cshear)))*VcRd;
          const vu=shearTorsionUtil(VplTRdi,Vx);
          if(vu>vt.u) vt={u:vu,x,V:Vx,T:Tt,tauT:tauTi,tauW:tauWi,VplTRd:VplTRdi,combo:se.combo.label,zeroCapacity:VplTRdi<=1e-9&&Vx>1e-9};
          tauT=Math.max(tauT,tauTi);
          tauW=Math.max(tauW,tauWi);
          // P385 3.1.2 cross-section interaction, Class 1/2 plastic, Class 3 elastic; M_z term = M_z,tot
          const u=cls12? Math.pow(My/Mply,2)+Mw/Mplf+MzTot/Mplz
                        : My/Mely+MzTot/Melz+Mw/Melf;
          return {x,My,phi,Mw,Mz,MzImp,MzTot,u};
        });
        rows.forEach(r2=>{ if(r2.u>cross.u) cross={...r2,combo:se.combo.label};
          MwMaxAbs=Math.max(MwMaxAbs,r2.Mw); MzMax=Math.max(MzMax,r2.Mz); MzTotMax=Math.max(MzTotMax,r2.MzTot); phiUmax=Math.max(phiUmax,Math.abs(r2.phi)); });
        g.p2.forEach((v,i)=>{ const Bi=Math.abs(EIw*v)/1e9; if(Bi>BMaxAbs){ BMaxAbs=Bi; BMaxPos=g.xs[i]; } });   // bimoment B = EI_w phi'', kN.m2
        grids.push({combo:se.combo,rows});
        const n=g.xs.length;
        const Tt0=O.GIt*g.p1[0]/1e6, TtL=O.GIt*g.p1[n-1]/1e6;    // kNm, St Venant part at the ends
        if(Math.abs(Tt0)>Math.abs(TtEnds[0])) TtEnds[0]=Tt0;
        if(Math.abs(TtL)>Math.abs(TtEnds[1])) TtEnds[1]=TtL;
        // total torque T = G I_T phi' - E I_w phi''' at the ends (St Venant + warping parts)
        const T0=(O.GIt*g.p1[0]-EIw*g.p3[0])/1e6, TL=(O.GIt*g.p1[n-1]-EIw*g.p3[n-1])/1e6;
        if(Math.abs(T0)>Math.abs(TEnds[0])) TEnds[0]=T0;
        if(Math.abs(TL)>Math.abs(TEnds[1])) TEnds[1]=TL;
      });
      const VplTRd=vt.VplTRd;                                                          // coincident V-T sweep, Eq 6.26/6.27
      // Method A comparison (simplified flange couple - conservative)
      let MwA=0; { let Tud=0;
        (O.sols[0]? O.sols[0].sol.xs:[]);
        // rebuild from the governing-moment combo torque set magnitudes
      }
      // simplified Mw: point torques -> flange SS BM; ud/lin -> F L/8
      const gm=O.sols.find(se=>se.combo===a.governM.combo)||O.sols[0];
      // Method of the elastic warping analysis (19 Sep 2026 gap closure, G4):
      // 'closed' = P385 App C Cases 3/4/10, 'fe' = the warping-torsion FE of
      // js/checks/torsion-fe.js (torsion cantilevers, partial-span torque,
      // warping-fixed ends). The FE mesh is doubled once; PASS is refused when
      // the change exceeds the tool threshold TORSION_FE_MESH_BLOCK.
      const feMethod=O.method==='fe';
      if(feMethod && !O.converged) unsupported.push('Warping-torsion FE mesh has not converged: doubling the mesh ('+O.nElemCoarse+' to '+O.nElem+' elements) changed the peak twist / St Venant torque / bimoment by '+(O.meshError*100).toFixed(2)+' % (limit '+(O.meshBlock*100).toFixed(1)+' %). The torsion effects are printed but PASS is blocked; refine the load layout or report the case.');
      // 20 Sep 2026 torsion + N/Mz: elastic yield criterion (6.1) with N_Ed, M_y, M_z,tot, sigma_w,
      // tau_V, tau_t, tau_w at the four section points (cl 6.2.7(5); helper torsionElasticOpen above)
      const elastic=torsionElasticOpen(a,sec,fy,gM0,O,F,MzImp);
      if(elastic && !(elastic.geom.Wn0>0 && elastic.geom.Sw1>0)) unsupported.push('Torsion: the P385 warping table (W_n0, S_w1'+(chan? ', W_n2, S_w2, S_w3' : '')+') is not tabulated for this section, so the warping stresses of the elastic yield criterion (6.1) cannot be evaluated; PASS is blocked.');
      tor={box:false,p385:true,TEd:a.tors.Tmax,governT:a.tors.governT,tp:sec.tp||null,
        elastic,MzImp,MzTot:MzTotMax,MzTwistMax:MzMax,
        method:O.method||'closed',methodLabel:O.methodLabel||'SCI P385 App C closed forms (Cases 3/4/10)',fe:feMethod,
        nElem:O.nElem||null,nElemCoarse:O.nElemCoarse||null,meshError:feMethod? O.meshError:null,meshBlock:O.meshBlock||null,meshConverged:feMethod? !!O.converged:true,
        bcText:O.bcText||'',feReasons:O.feReasons||[],
        BMax:BMaxAbs,BMaxPos,TEnds,nSolves:(O.nSolves!=null? O.nSolves : null),nCached:(O.nCached!=null? O.nCached : null),
        e0:(sec.tp&&sec.tp.e0!=null)? sec.tp.e0 : (sec.e0!=null? sec.e0*10 : null),
        esc:(sec.tp&&sec.tp.esc!=null)? sec.tp.esc : null,
        aa:O.aa,X:O.X,IT:O.IT,Iw:O.Iw,chan,cls12,
        Mply,Mplz,Mplf,Mely,Melz,Melf,
        cross,grids,MwMax:MwMaxAbs,MzMax,phiUmax,
        TtEnds,tauT,tauW,VplTRd,vt,vtUtil:vt.u, vtZeroCapacity:!!vt.zeroCapacity,
        phiSer:O.sls? Math.abs(O.sls.phiMax):0, phiSerDeg:O.sls? Math.abs(O.sls.phiMax)*180/Math.PI:0,
        phiSerPos:O.sls? O.sls.phiPos/1000:0, governTw:O.sls? O.sls.combo.label:'',
        GIt:O.GIt};
    } else {
      unsupported.push("Torsion on this open section is NOT COVERED: "+((a.torsO&&a.torsO.reason)||"the P385 St Venant + warping analysis requires a warping constant")+".");
      tor={box:false,TEd:a.tors.Tmax,governT:a.tors.governT,tp:sec.tp||null,
        e0:(sec.tp&&sec.tp.e0!=null)? sec.tp.e0 : (sec.e0!=null? sec.e0*10 : null),
        esc:(sec.tp&&sec.tp.esc!=null)? sec.tp.esc : null};
    }
  }
  if(a.torsErr) unsupported.push("Eccentric loads are active but torsion cannot be evaluated: "+a.torsErr+".");
  // 20 Sep 2026 torsion + N/Mz: the former block "Combined torsion with direct axial force or
  // imposed minor-axis bending is not implemented as one interaction" is removed; the
  // verification is the (6.1) elastic check above, 6.61/6.62 with M_z,Ed = M_z + phi.M_y
  // (annexB2) and EN 1993-6 (A.1) with M_z,tot (unrestrained paths); the basis is printed
  // (torsionCombinedBasis, set per path once buck / annex are known - below and in the
  // unrestrained paths). 20 Sep 2026 review: the binding policy of (6.1) (elasticBindingPolicy).
  if(tor && (tor.p385 || tor.box)){
    tor.combinedActive=(Math.abs(F)>1e-9 || MzEd>1e-9);
    if(!tor.elastic) unsupported.push('Torsion: the elastic yield criterion (6.1) could not be evaluated (no coincident station data); PASS is blocked.');
    else { const pol=elasticBindingPolicy(sec,cl,F); tor.elastic.binding=pol.binding; tor.elastic.bindingBasis=pol.basis; }
  }
  // SLS twist guideline: SCI P385 suggests limiting the serviceability rotation
  // to about 2 degrees; flagged as a non-blocking advisory (guideline, not a
  // code limit) - matches common commercial-software practice.
  {
    const twistDeg = tor? (tor.p385? tor.phiSerDeg : (tor.box? tor.phiDeg : null)) : null;
    if(twistDeg!=null && twistDeg>2)
      advisory.push("Serviceability twist &theta;<sub>max</sub> = "+twistDeg.toFixed(1)+"&deg; exceeds the ~2&deg; guideline (SCI P385): check that the rotation is acceptable for whatever the member supports.");
  }
  // moment resistance, cl 6.2.5, with the cl 6.2.8 shear check at the point of maximum moment
  const Zx=sec.Zx*1e3, Sx=sec.Sx*1e3;
  const Wy = cl.cls<=2? Sx : Zx;
  let McRd=Wy*fy/gM0/1e6;                                     // kN.m
  const gfb=a.governM.fb;
  // kN, worst shear coincident with a combination's own maximum moment - swept
  // over EVERY enabled ULS combination, not just the governing-moment one (a
  // combo with a slightly lower Mmax can pair it with a much higher V there).
  const VatM=Math.max(...a.ulsResults.map(res=>Math.abs(interpAt(res.fb.xs,res.fb.V,res.Mpos))))/1000;
  const VplMoment=(tor&&tor.VplTRd!=null)? tor.VplTRd : VcRd;
  const halfVpl=0.5*VplMoment; // cl 6.2.8(4): use Vpl,T,Rd when torsion is present
  const lowShearAtM = VatM<=halfVpl+1e-9;
  let hsNote=null, mvForm=null, rhoAtM=null;
  if(!lowShearAtM){
    if(VatM>VplMoment+1e-9){
      hsNote="V<sub>Ed</sub> at the point of maximum moment exceeds V<sub>pl,Rd</sub>: the member has already failed the pure shear resistance check, so the cl 6.2.8 reduced moment formula is not valid";
    } else if(cl.cls<=3){
      // cl 6.2.8(3): (1 - rho) f_y on the shear area, every family and class
      // (shearReducedResistances): rolled I/H Eq 6.30 form, elastic I/H Class 3,
      // channel and hollow-section web moduli
      rhoAtM=Math.min(Math.pow(2*VatM/VplMoment-1,2),1);
      const R=shearReducedResistances(sec,cl,fy,Av,rhoAtM,gM0);
      McRd=Math.min(McRd,R.MvY); mvForm=R.form;
      hsNote="V<sub>Ed</sub> at the point of maximum moment exceeds 0.5V<sub>pl,Rd</sub>: moment resistance reduced per cl 6.2.8(3) with &rho;=(2V<sub>Ed</sub>/V<sub>pl,Rd</sub>&minus;1)&sup2; = "+g(rhoAtM,3)+", M<sub>v,Rd</sub> = "+R.form;
    }
  }
  const Mx=Math.abs(a.Mmax);
  const momUtil=McRd>0? Mx/McRd : 0;
  // ---- span-wise coexistent M-V (cl 6.2.8) and M-V-N / M-V-Mz (cl 6.2.10) sweep ----
  // Every station of every enabled ULS combination's own (V, M) fields (the
  // interaction is nonlinear, so the envelopes of V and M do not bound it).
  // Where V > 0.5 V_pl(,T),Rd: rho = (2V/V_pl - 1)^2 and the resistances of
  // shearReducedResistances() apply at that station - all families, Class 1-3:
  //   coex : M/M_v,Rd (6.2.8), worst station;
  //   mvn  : with N_Ed != 0 or M_z != 0 (6.2.10(3)): the 6.2.9 interaction
  //          re-evaluated with N_V,Rd, M_v,y,Rd, M_v,z,Rd and a_V -
  //          I/H and RHS Class 1/2: M_N,V,y,Rd = M_v,y,Rd (1 - n_V)/(1 - 0.5 a_V)
  //          <= M_v,y,Rd with the 6.2.9.1(4) waiver (N <= 0.25 N_V,Rd and N <=
  //          0.5 h_w t_w (1 - rho) f_y), M_N,V,z,Rd per 6.2.9.1(5) / Eq 6.40;
  //          uniaxial: M_y/M_N,V,y,Rd; biaxial: (M_y/M_N,V,y)^alpha + (M_z/M_N,V,z)^beta;
  //          Class 3 and channels: n_V + M_y/M_v,y,Rd + M_z/M_v,z,Rd (6.2.9.2 / 6.2.1(7)).
  // V > V_pl,Rd at a station is a pure shear failure (the reduced-moment forms
  // are not valid there) and is reported as such.
  let coex=null, mvn=null;
  if(cl.cls<=3){
    const VplB=(tor&&tor.VplTRd!=null)? tor.VplTRd : VcRd;
    const NEdAbs=Math.abs(F), hasN=NEdAbs>1e-9, hasMz=MzEd>1e-9;
    const hw62=sec.D-2*sec.tf;
    let worst={u:momUtil,x:a.Mpos*1000,red:false};
    let worstN={u:-1};
    a.ulsResults.forEach(res=>{ const cfb=res.fb;
    cfb.xs.forEach((x,i)=>{
      const Vx=Math.abs(cfb.V[i])/1e3, Mxx=Math.abs(cfb.M[i])/1e6;
      if(Mxx<1e-9 && !hasN && !hasMz) return;
      if(Vx>VplB+1e-9){
        const u=Vx/Math.max(VplB,1e-9);
        if(!worst.pureShearFail || Mxx>worst.M+1e-9 || (Math.abs(Mxx-worst.M)<=1e-9 && u>worst.u))
          worst={u,x,red:true,pureShearFail:true,V:Vx,M:Mxx,VplRd:VplB,MvRd:null,combo:res.combo.label};
        return;
      }
      if(!(Vx>0.5*VplB)) return;          // low shear: M/M_c,Rd <= momUtil already (M_c,Rd carries the peak-station reduction)
      const rho=Math.min(Math.pow(2*Vx/VplB-1,2),1);
      const R=shearReducedResistances(sec,cl,fy,Av,rho,gM0);
      const MvRd=Math.min(McRd,R.MvY);
      if(Mxx>1e-9){
        const u=Mxx/Math.max(MvRd,1e-9);
        if(!worst.pureShearFail && u>worst.u) worst={u,x,red:true,V:Vx,M:Mxx,MvRd,rho,form:R.form,combo:res.combo.label};
      }
      if(hasN||hasMz){
        const nV=NEdAbs/Math.max(R.NV,1e-9);
        let MNVy, MNVz, alpha=1, beta=1, form, waiver=false, aV=R.aV;
        if(cl.cls<=2 && sec.kind==='I'){
          MNVy=Math.max(0,Math.min(R.MvY*(1-nV)/(1-0.5*aV),R.MvY));
          if(NEdAbs<=0.25*R.NV && NEdAbs*1000<=0.5*hw62*sec.tw*(1-rho)*fy/gM0){ MNVy=R.MvY; waiver=true; }
          alpha=2; beta=Math.max(5*nV,1);
          MNVz= nV<=aV? R.MvZ : R.MvZ*(1-Math.pow((nV-aV)/(1-aV),2));
          form='plastic, 6.2.9.1 with (1 &minus; &rho;)f<sub>y</sub> on A<sub>v</sub>';
        } else if(cl.cls<=2 && sec.isBox){
          const afV=R.afV;
          MNVy=Math.max(0,Math.min(R.MvY*(1-nV)/(1-0.5*aV),R.MvY));
          MNVz=Math.max(0,Math.min(R.MvZ*(1-nV)/(1-0.5*afV),R.MvZ));
          alpha= nV<=0.8? 1.66/(1-1.13*nV*nV) : 6; beta=alpha;
          form='plastic, Eq 6.39/6.40 with (1 &minus; &rho;)f<sub>y</sub> on A<sub>v</sub>';
        } else {
          MNVy=R.MvY; MNVz=R.MvZ;
          form= sec.kind==='channel'? 'linear, cl 6.2.1(7) with N<sub>V,Rd</sub>, M<sub>v,y,Rd</sub>, M<sub>v,z,Rd</sub>' : 'elastic, cl 6.2.9.2 with N<sub>V,Rd</sub>, M<sub>v,y,Rd</sub>, M<sub>v,z,Rd</sub>';
        }
        let u;
        if(cl.cls<=2 && sec.kind!=='channel'){
          u= hasMz? Math.pow(Mxx/Math.max(MNVy,1e-9),alpha)+Math.pow(MzEd/Math.max(MNVz,1e-9),beta) : Mxx/Math.max(MNVy,1e-9);
        } else {
          u= nV + Mxx/Math.max(MNVy,1e-9) + (hasMz? MzEd/Math.max(MNVz,1e-9) : 0);
        }
        if(u>worstN.u) worstN={u,x,V:Vx,M:Mxx,Mz:MzEd,N:NEdAbs,rho,NV:R.NV,nV,MvY:R.MvY,MvZ:R.MvZ,MNVy,MNVz,alpha,beta,aV,afV:R.afV,waiver,form,formY:R.form,formZ:R.formZ,combo:res.combo.label,
          biax:hasMz,plastic:(cl.cls<=2 && sec.kind!=='channel')};
      }
    });
    });
    if(worst.red){
      coex=worst;
      if(!hsNote) hsNote=worst.pureShearFail
        ? "Coexistent shear and moment (cl 6.2.8): at x = "+(worst.x/1000).toFixed(2)+" m, V<sub>Ed</sub> = "+worst.V.toFixed(0)+" kN exceeds V<sub>pl,Rd</sub> = "+worst.VplRd.toFixed(0)+" kN, so pure shear failure governs and the reduced M<sub>v,Rd</sub> formula is not valid"
        : "Coexistent shear and moment (cl 6.2.8): at x = "+(worst.x/1000).toFixed(2)+" m, V<sub>Ed</sub> = "+worst.V.toFixed(0)+" kN &gt; 0.5V<sub>pl,Rd</sub> and the reduced M<sub>v,Rd</sub> = "+worst.MvRd.toFixed(0)+" kN&middot;m governs";
    }
    if(worstN.u>=0 && !(coex&&coex.pureShearFail)) mvn=worstN;
  }
  // vertical deflection (NA 2.23) - governing enabled SLS combination and segment;
  // the limit is the segment's own (span/divisor, L/divisorCant for a cantilever
  // segment, capped by the optional absolute limit) as set by analyse()
  const span=a.deflection?a.deflection.span:a.L;
  const divisor=(a.deflection&&a.deflection.divisor!=null)? a.deflection.divisor : S.divisor;
  const dlimit=(a.deflection&&a.deflection.limit!=null)? a.deflection.limit : span/divisor;
  const dmax=Math.abs(a.deflection?a.deflection.dmax:a.dmax), defOk=dmax<=dlimit;
  const deflCant=!!(a.deflection&&a.deflection.cant), deflAbsGoverns=!!(a.deflection&&a.deflection.absGoverns);
  // uplift / hold-down at every support, every combination (EN 1990 2.4.4)
  const holdDown=holdDownCheck(a);
  holdDown.unsupported.forEach(m=>unsupported.push(m));
  holdDown.advisory.forEach(m=>advisory.push(m));
  if(a.stability&&a.stability.notes) a.stability.notes.forEach(m=>advisory.push(m));
  if(a.companionNote) advisory.push(a.companionNote);
  const isCantR=isCantilever(S);
  const buck=(ax && !ax.tension)? annexB2(a,sec,fy,cl,McRd,true,isCantR,aeff) : null; // fully restrained: not susceptible -> Table B.1; MbRd = Mc,Rd
  if(buck && buck.tfb && !buck.tfb.ok) unsupported.push('PFC under axial compression: '+buck.tfb.reason+'; torsional / torsional-flexural buckling (cl 6.3.1.4) cannot be verified, PASS is blocked.');
  const utils=[
    {name:"Shear  V_Ed/V_c,Rd",val:shearUtil},
    {name:"Bending  M_Ed/M_c,Rd",val:momUtil},
    {name:"Deflection",val:dmax/dlimit},
  ];
  if(ax){
    utils.push({name: ax.tension? "Tension  N_Ed/N_t,Rd" : (ax.aeff&&ax.aeff.active? "Compression  N_Ed/N_c,Rd (A_eff)" : "Compression  N_Ed/N_pl,Rd"), val:ax.nUtil});
    utils.push({name: ax.biax? ("Biaxial bending"+(Math.abs(F)>1e-9?" + axial":"")+" (6.2.9.1)") : "Bending+axial cross-section (6.2.9)",val:ax.mUtil});
    if(!ax.tension && buck && buck.Fc>1e-9){   // cl 6.3.3 applies only with axial compression
      if(buck.tfb && buck.tfb.ok) utils.push({name:TFB_UTIL_NAME,val:buck.tfb.util});
      utils.push({name:"Member buckling y-y (Eq 6.61)",val:buck.u1});
      utils.push({name:"Member buckling z-z (Eq 6.62)",val:buck.u2});
    }
  }
  if(coex) utils.push({name:coex.pureShearFail? "Pure shear failure at M-V check point (6.2.6)" : "Bending+shear coexistent (6.2.8)",val:coex.u});
  if(mvn) utils.push({name:mvnUtilName(mvn),val:mvn.u});
  if(web&&web.checked){
    utils.push({name:WEB_UTIL_NAMES[0],val:web.util2});
    utils.push({name:WEB_UTIL_NAMES[1],val:web.util72});
  }
  if(tor&&tor.box){
    utils.push({name:"Torsion  T_Ed/T_Rd",val:tor.torUtil});
    utils.push({name:"Shear+torsion  V_Ed/V_pl,T,Rd",val:tor.vtUtil});
  }
  if(tor&&tor.p385){
    utils.push({name:"Bending+torsion cross-section (P385 3.1.2)",val:tor.cross.u});
    utils.push({name:"Shear+torsion  V_Ed/V_pl,T,Rd",val:tor.vtUtil});
  }
  // 20 Sep 2026 torsion + N/Mz: elastic yield criterion (6.1), cl 6.2.7(5), computed in every torsion case;
  // 20 Sep 2026 review: verdict-binding (utils) only where the code gives no plastic route (Class 3, or an
  // open Class 1/2 section with N_Ed - elasticBindingPolicy), otherwise information (c.info; an advisory
  // when it exceeds 1 - 6.2.7(5) is permissive and 6.2.7(6) admits the plastic route)
  const info=[];
  if(tor&&tor.elastic){
    if(tor.elastic.binding) utils.push({name:ELASTIC_TORSION_UTIL_NAME,val:tor.elastic.u});
    else {
      info.push({name:ELASTIC_TORSION_UTIL_NAME,val:tor.elastic.u,note:tor.elastic.bindingBasis});
      if(tor.elastic.u>1.0001) advisory.push('ADVISORY - elastic yield criterion (6.1) with torsion, cl 6.2.7(5) = '+tor.elastic.u.toFixed(3)+' at x = '+(tor.elastic.x/1000).toFixed(2)+' m ('+tor.elastic.point+'): first yield is reached under the design actions (information only, not a utilisation: '+tor.elastic.bindingBasis+').');
    }
  }
  // 20 Sep 2026 torsion + N/Mz: advisory superposition of Eq 6.62 and (A.1) (fully restrained: k_alpha = 1);
  // the unrestrained paths recompute it with their own k_alpha and print their own text
  if(tor&&tor.p385){ tor.superposition=torsionSuperposition(tor,null,buck); if(tor.superposition && (S.restraint||'full')==='full') advisory.push(tor.superposition.text); }
  // 20 Sep 2026 review: the basis text of this (fully restrained) path; the unrestrained paths overwrite it with theirs
  if(tor&&(tor.p385||tor.box)) tor.combinedBasis=torsionCombinedBasis({box:!!tor.box,tension:!!(ax&&ax.tension),hasN:Math.abs(F)>1e-9,hasMz:MzEd>1e-9,restrained:true,buckEvaluated:!!(buck&&buck.Fc>1e-9),annexEvaluated:false,binding:!!(tor.elastic&&tor.elastic.binding)});
  let gov=utils[0]; utils.forEach(u=>{ if(u.val>gov.val) gov=u; });
  const pass=unsupported.length===0 && utils.every(u=>u.val<=1.0001);
  return {sci:true,tor,coex,mvn,web,ax,aeff,buck,eps,cl,clsName,unsupported,advisory,info,fy,eta,hw,cOut,Av,AvRaw,avFloor,VcRd,Fv,shearUtil,
    sbRatio,sbLimit,sbOk,Zx,Sx,Wy,McRd,hsNote,mvForm,rhoAtM,VatM,halfVpl,lowShearAtM,Mx,momUtil,F,
    span,divisor,dlimit,dmax,defOk,deflCant,deflAbsGoverns,holdDown,restraintForces:null,utils,gov,pass};
}

/* ===========================================================================
   STANDARD (closed-form) Mcr METHOD  -  S.mcrMethod === 'standard'
   ---------------------------------------------------------------------------
   Everything from here to the end of checksEC3UnrestrainedSCI() is the
   standard-method implementation: C1 from the NCCI SN003a tables / SCI
   end-moment curve / Serna quarter-point expression, the SN003a closed-form
   Mcr with the C2*zg load-height term where C2 is published, the SN006a
   cantilever C factors (sn006C in 01-computation-engine.js), and the P385/P362
   channel kappa chain. This is how reference software-type software derives Mcr.
   js/08-mcr-eigen-patch.js keeps this function alive as
   window.checksEC3UnrestrainedStandard and delegates to it when the user
   selects the standard method; it is NOT dead code. The helpers it depends on
   (C1_END_MOMENT, SN006/sn006C, sernaC1, c1FromPsi, computeC1, mcrEC3,
   cmTableB3, annexB2) must stay in place.
   =========================================================================== */

// ---- Moment ordinate at an interior station of a diagram with jumps ----
// An applied in-span couple makes M(x) discontinuous. A quarter-point or
// mid-span sample that lands exactly on the jump must take the larger of the
// two side ordinates: the closed-form C1 expressions (Serna quarter points,
// the reference software M_o) describe the envelope of the diagram, and the far-side
// value overstates C1 (19 Sep 2026 verification campaign, UB-27: C1 1.68 ->
// 1.21, the standard M_cr was 20 % above the eigenvalue). The diagram grid
// carries x - 1e-4, x, x + 1e-4 at a jump, so 1e-3 mm either side reaches
// both sides while a smooth diagram is unchanged to 1e-6. Clamped to the
// segment. Pure.
function mAtStation(fb,x,xa,xb){
  const d=1e-3;
  const l=interpAt(fb.xs,fb.M,Math.max(x-d,xa)), r=interpAt(fb.xs,fb.M,Math.min(x+d,xb));
  return Math.abs(l)>=Math.abs(r)? l : r;
}
// ---- C1 inputs in the reference software convention (pure) ----
// fb = {xs (mm), M (N.mm)} of one combination's BMD (sagging positive);
// xa, xb = the segment ends (mm). Returns kN.m:
//   M1, M2 = the end moments of the segment, M2 the larger in magnitude;
//   Mo     = mid-segment moment above the chord joining M1 and M2 (the free
//            bending moment from the loads inside the segment);
//   psi    = M1/M2 (clamped to [-1, 1]);
//   mu     = Mo/M2, capped at +/-300 (reference software prints 300.000 when M2 ~ 0).
// End moments are read a fraction inside the segment so that a point moment or
// a support reaction exactly at the end does not pick the wrong side of the jump.
function c1Inputs(fb,xa,xb){
  const Ma=interpAt(fb.xs,fb.M,xa+1e-4)/1e6, Mb=interpAt(fb.xs,fb.M,xb-1e-4)/1e6;
  const Mmid=mAtStation(fb,(xa+xb)/2,xa,xb)/1e6;
  let [M1,M2]= Math.abs(Mb)>=Math.abs(Ma)? [Ma,Mb] : [Mb,Ma];
  const Mo=Mmid-(Ma+Mb)/2;
  // end moments that are numerical noise next to the in-span moment (a pinned
  // end) are reported as zero, so psi and mu are deterministic
  const scale=Math.max(Math.abs(M1),Math.abs(M2),Math.abs(Mo),1e-9);
  if(Math.abs(M2)<=1e-6*scale){ M1=0; M2=0; }
  else if(Math.abs(M1)<=1e-6*scale) M1=0;
  const eps=1e-9;
  const psi= Math.abs(M2)>eps? Math.max(-1,Math.min(1,M1/M2)) : 1;
  let mu= Math.abs(M2)>eps? Mo/M2 : (Math.abs(Mo)>eps? 300*Math.sign(Mo) : 0);
  mu=Math.max(-300,Math.min(300,mu));
  return {M1,M2,Mo,psi,mu,xa,xb};
}
// ---- Governing segment for the C1 inputs ----
// Whole member for a span without intermediate lateral restraints or a
// cantilever. With intermediate restraints: the bay between adjacent
// lateral-restraint points (ends holding U_y and v-restraints) that contains
// the governing combination's peak moment. Returns {xa, xb, whole} in mm.
function c1Segment(a){
  const L=a.L;
  const isCant=isCantilever(S);
  const pts=[...new Set(endsList().filter(e=>e.uy).map(e=>+((e.x*1000).toFixed(3)))
    .concat((S.ltbRestraints||[]).filter(r=>r.v!==false).map(r=>+(((+r.pos)*1000).toFixed(3))))
    .filter(x=>isFinite(x)&&x>=-1e-6&&x<=L+1e-6))].sort((p,q)=>p-q);
  if(isCant || pts.length<3) return {xa:0,xb:L,whole:true};
  const xm=a.Mpos*1000;
  for(let i=0;i<pts.length-1;i++) if(xm>=pts[i]-1e-6 && xm<=pts[i+1]+1e-6) return {xa:pts[i],xb:pts[i+1],whole:false};
  if(xm<pts[0]) return {xa:0,xb:pts[0],whole:false};
  return {xa:pts[pts.length-1],xb:L,whole:false};
}
// ---- LTB buckling curve, UK NA Table NA.1 (cl 6.3.2.3), shared by both Mcr methods ----
// Rolled doubly symmetric I/H sections and hot-finished hollow sections:
// h/b <= 2 -> b, 2 < h/b <= 3.1 -> c, h/b > 3.1 -> d. Cold-formed hollow
// sections (with welded sections): h/b <= 2 -> c, otherwise d. Channels (not
// doubly symmetric, "all other hot-rolled sections") -> d.
function ltbCurveNA(sec){
  if(sec.isBox && sec.boxType==='CF') return sec.D/sec.B<=2? {alphaLT:0.49,curve:'c'} : {alphaLT:0.76,curve:'d'};
  if(sec.kind==='channel') return {alphaLT:0.76,curve:'d'};
  const hb=sec.D/sec.B;
  return hb<=2? {alphaLT:0.34,curve:'b'} : hb<=3.1? {alphaLT:0.49,curve:'c'} : {alphaLT:0.76,curve:'d'};
}
// ---- SN003a closed-form Mcr (pure), N.mm ----
// Doubly symmetric section, k = kw = 1, G = 81000 N/mm2 (SN003a / P385):
// Mcr = C1 (pi^2 E Iz/LE^2) { sqrt[ Iw/Iz + LE^2 G It/(pi^2 E Iz) + (C2 zg)^2 ] - C2 zg }
// The C2 zg term is applied only when zg is non-zero AND C2 is published for
// the recognised diagram (C2 = null otherwise). Hollow sections have Iw = 0.
function mcrClosedForm(sec,E,LE,C1,C2,zg){
  const G=81000, Iz=sec.Iy*1e4, It=sec.J*1e4, Iw=(sec.Iw||0)*1e12;
  const T1=Math.PI*Math.PI*E*Iz/(LE*LE);          // N
  const IwIz=Iw/Iz;                                // mm2
  const GIt=G*It;                                  // N.mm2
  const zgUsed=(Math.abs(zg||0)>1e-9 && C2!=null && C2>0);
  const zgTerm=zgUsed? C2*zg : 0;
  const Mcr=C1*T1*(Math.sqrt(Math.max(IwIz+GIt/T1+zgTerm*zgTerm,0))-zgTerm);
  return {Mcr,T1,IwIz,GIt,zgUsed,zgTerm};
}
// ---- Load height for the closed-form route (pure, reads S) ----
// The SN003a form takes ONE z_g for the whole diagram. It is the most
// destabilising signed height among the transverse loads that are active in
// the combination `factors`: each load's own z_g from loadZgValue() (the
// per-load value when eccOn is on, else the default S.za), positive above
// the shear centre for a DOWNWARD load. A load acting upward at +z_g is
// stabilising, so its height enters with the sign reversed. Automatic
// self-weight (centroid, z_g = 0) and moment loads carry no load height.
// Falls back to S.za when no transverse load is active in the combination.
// `combo` is a combination object ({factors, mask?}) or a bare factors object;
// pattern combinations see only the loads their mask keeps (comboLoadPieces).
// Returns {zg (mm), source: 'loads' | 'default'}.
function asCombo(combo){ return (combo && combo.factors)? combo : {factors:combo||{}}; }
function stdZgFor(combo){
  let best=null;
  comboLoadPieces(asCombo(combo)).forEach(p=>{
    if(p.type==='moment') return;
    const f=p.factor; if(!f) return;
    const mag= p.type==='point'? p.P : (p.w1+p.w2)/2;
    if(Math.abs(mag)<1e-12) return;
    const z0=(typeof loadZgValue==='function')? loadZgValue(p.ld) : (+S.za||0);
    const zg=(mag*f<0)? -z0 : z0;
    if(best==null||zg>best) best=zg;
  });
  return best==null? {zg:(+S.za||0), source:'default'} : {zg:best, source:'loads'};
}
// ---- Transverse-load shape of one combination (pure, reads S) ----
// Counts the loads with a non-zero factor (self-weight excluded, it is always
// a full-span UDL): full-span UDLs (trap with w1 = w2 counts), point loads
// (and how many of them sit at mid-span, +/-1% of L), any other transverse
// load, and applied moments. The SN003a Table 3.2 rows are recognised from
// these counts, not from the quarter-point moment ratio alone.
function stdLoadShape(combo){
  const L=(+S.L)*1000;
  let nUdl=0,nPoint=0,nCentral=0,nOther=0,nMoment=0;
  comboLoadPieces(asCombo(combo)).forEach(p=>{
    const f=p.factor; if(!f) return;
    if(p.type==='moment'){ if(Math.abs(p.M)>1e-12) nMoment++; return; }
    if(p.type==='point'){
      if(Math.abs(p.P)<1e-12) return;
      nPoint++; if(Math.abs(p.pos-L/2)<=0.01*L) nCentral++; return;
    }
    const w1=p.w1, w2=p.w2;
    if(Math.abs(w1)<1e-12&&Math.abs(w2)<1e-12) return;
    const full=p.x1<=1e-3 && Math.abs(p.x2-L)<=1e-3 && Math.abs(w1-w2)<=1e-9*Math.max(Math.abs(w1),1);
    if(full) nUdl++; else nOther++;
  });
  return {nUdl,nPoint,nCentral,nOther,nMoment,
    notLoaded: nUdl+nPoint+nOther===0,                                   // only applied moments (and self-weight)
    uniform: nMoment===0 && nUdl>0 && nPoint===0 && nOther===0,           // full-span UDL family only
    central: nMoment===0 && nPoint>0 && nCentral===nPoint && nUdl===0 && nOther===0}; // central point load(s) only
}
// ---- Status of the entered load height on the closed-form route (pure, reads S.destab) ----
// Returns {text (HTML, always printed), block (HTML or null: blocks PASS),
// advisory (HTML or null)}. The SN003a C2 z_g term exists only for the four
// Table 3.2 rows; a destabilising height on any other diagram cannot be
// included by the closed form, and silently taking the load at the shear
// centre would be unconservative, so PASS is blocked unless the user has
// chosen the destabilising L_E = 1.2 L treatment (BS 5950 practice) instead.
function stdZgStatus(zg,C2,zgUsed){
  zg=+zg||0;
  const zs=(zg>0?'+':'')+zg.toFixed(0)+' mm';
  if(Math.abs(zg)<1e-9) return {text:'z<sub>g</sub> = 0 (load through the shear centre)', block:null, advisory:null};
  if(zgUsed){
    return {text:'z<sub>g</sub> = '+zs+' above the shear centre, C<sub>2</sub> = '+(+C2).toFixed(3)+' (SN003a Table 3.2)'+(S.destab? '; the destabilising &times;1.2 L<sub>E</sub> switch is also on' : ''),
      block:null,
      advisory: S.destab? 'Standard (closed-form) M<sub>cr</sub>: the destabilising &times;1.2 L<sub>E</sub> switch AND the SN003a C<sub>2</sub>z<sub>g</sub> load-height term are both applied. They describe the same effect, so the load height is counted twice (conservative). Untick the switch or set z<sub>g</sub> = 0 to apply one treatment only.' : null};
  }
  if(zg>0 && S.destab) return {text:'z<sub>g</sub> = '+zs+' entered &mdash; NOT applied in M<sub>cr</sub> (C<sub>2</sub> not published for this moment diagram); the load height is carried by the destabilising L<sub>E</sub> = 1.2&times;L<sub>E</sub>-factor&times;L only (BS 5950 practice)',
    block:null,
    advisory:'Standard (closed-form) M<sub>cr</sub>: z<sub>g</sub> = '+zs+' (destabilising) is entered but SN003a publishes C<sub>2</sub> only for the simply supported and fixed-ended full-span UDL and central point-load diagrams. The C<sub>2</sub>z<sub>g</sub> term is therefore NOT applied; the load height is represented by the destabilising L<sub>E</sub> = 1.2&times;L<sub>E</sub>-factor&times;L switch only. Use the FE eigenvalue method for the exact load-height effect.'};
  if(zg>0) return {text:'z<sub>g</sub> = '+zs+' entered &mdash; NOT applied: C<sub>2</sub> not published for this moment diagram; M<sub>cr</sub> is at the shear centre (PASS blocked)',
    block:'Standard (closed-form) M<sub>cr</sub>: a destabilising load height z<sub>g</sub> = '+zs+' is entered, but SN003a publishes C<sub>2</sub> only for the simply supported and fixed-ended full-span UDL and central point-load diagrams (Table 3.2). The closed form cannot include the load height for this moment diagram, so the printed M<sub>cr</sub> is computed with the load at the shear centre (unconservative). Use the FE eigenvalue M<sub>cr</sub> method, a verified M<sub>cr</sub>, or the destabilising L<sub>E</sub> switch (BS 5950 practice); PASS is blocked.',
    advisory:null};
  return {text:'z<sub>g</sub> = '+zs+' (below the shear centre, stabilising) &mdash; not applied: C<sub>2</sub> not published for this moment diagram; M<sub>cr</sub> is taken with the load at the shear centre (conservative)', block:null, advisory:null};
}
function sn003aC1(a,isCant,seg,diag){
  // C1 for the Mcr calculation (SN003a Table 3.1/3.2) and 1/sqrt(C1) for the
  // P362 Eq 6.55 simplified slenderness / NA 2.18 kc factor.
  // `seg` = {xa, xb} (mm) evaluates the diagram over that segment; omitted =
  // the whole member, using the analysis quantities exactly as before.
  // `diag` = {fb, factors}: the moment diagram the C1 describes and the load
  // factors of the combination that produced it (used for the load-list shape
  // recognition); omitted = the governing-moment combination a.governM.
  // Returns {C1, C2, label, route, c1in}; route is one of 'override',
  // 'cantilever', 'negligible', 'end-moment', 'uniform', 'point',
  // 'fixed-uniform', 'fixed-point', 'serna'.
  const fbM=(diag&&diag.fb)||a.governM.fb;
  const fac=(diag&&diag.combo)||((diag&&diag.factors)? {factors:diag.factors} : a.governM.combo);   // combination (pattern-aware)
  const xa=seg? seg.xa : 0, xb=seg? seg.xb : a.L, Ls=xb-xa;
  const c1in=c1Inputs(fbM,xa,xb);
  if(S.C1o!=null) return {C1:S.C1o, C2:null, label:'user override', route:'override', c1in};
  if(isCant) return {C1:1.0, C2:null, label:'cantilever', route:'cantilever', c1in};
  const Mat=x=>interpAt(fbM.xs,fbM.M,x)/1e6;
  const Mst=x=>mAtStation(fbM,x,xa,xb)/1e6;   // interior stations: larger side of a jump (in-span couple)
  // peak, end (read a fraction inside the segment) and quarter-point moments of
  // the segment; for the whole member these equal the analysis quantities
  let Mm=0; fbM.xs.forEach((x,i)=>{ if(x>=xa-1e-6&&x<=xb+1e-6) Mm=Math.max(Mm,Math.abs(fbM.M[i])/1e6); });
  const M0=Mat(xa+1e-4), ML=Mat(xb-1e-4), Mq=Mst(xa+Ls/4), Mq3=Mst(xa+3*Ls/4);
  if(Mm<1e-9) return {C1:1.0, C2:0, label:'negligible moment', route:'negligible', c1in};
  const endLevel=Math.max(Math.abs(M0),Math.abs(ML))/Mm;
  const shape=stdLoadShape(fac);
  // The end-moment curve ("Not Loaded") applies when the combination carries
  // no transverse load other than self-weight (decided from the load list, so
  // a light section's factored self-weight curvature does not disqualify it)
  // or when the diagram IS a straight line between the segment ends (every
  // grid value within 5% of Mmax of the chord); in both cases the larger end
  // moment must be the peak of the segment.
  const isLinear=fbM.xs.every((x,i)=> x<xa+1e-4 || x>xb-1e-4 || Math.abs(fbM.M[i]/1e6-(M0+(ML-M0)*(x-xa)/Ls))<=0.05*Mm);
  if(endLevel>0.98 && (shape.notLoaded || isLinear)){
    // psi = smaller end moment / larger end moment (SN003a Table 3.1 convention;
    // reference software psi = M1/M2), so a larger moment at x = L does not clamp to 1.
    const psi=c1in.psi;
    const c=Math.pow(1.33-0.33*psi,2); // SCI curve C1=(1.33-0.33psi)^2 = 1.77-0.88psi+0.11psi^2 (NA kc inverted)
    return {C1:c, C2:0, label:'linear end-moment gradient, &psi; = '+psi.toFixed(2)+' (SCI curve, NA 2.18)'+(shape.notLoaded? '; not loaded' : ''), route:'end-moment', c1in};
  }
  // The four tabulated transverse-load shapes (SN003a Table 3.2, k = kw = 1)
  // are recognised from the LOAD LIST of the combination (full-span UDL family
  // only, or central point load(s) only, no applied moments) on the whole
  // member with a vertical support (U_z) at each end and no hinge inside it:
  // simply supported ends (R_y free) with no significant end moment, or both
  // ends fixed (U_z + R_y). The quarter-point/mid-span ratio of the diagram
  // is kept as a guard (the point-load row tolerates a self-weight moment
  // share of about 9%).
  const [eA,eB]=endsList();
  const endSupported=!!(eA.uz&&eB.uz);
  const endFixed=!!(eA.uz&&eB.uz&&eA.ry&&eB.ry);
  const interiorBreak=(S.hinges||[]).some(h=>(+h.pos)*1000>xa+1e-6&&(+h.pos)*1000<xb-1e-6);
  const whole=(!seg)||(seg.xa<=1e-6&&Math.abs(seg.xb-a.L)<=1e-6);
  const r=(Math.abs(Mq)+Math.abs(Mq3))/(2*Mm);
  if(whole && endSupported && !interiorBreak){
    if(!endFixed && endLevel<0.02){
      if(shape.uniform && Math.abs(r-0.75)<=0.02) return {C1:1.127, C2:0.454, label:'simply supported + uniformly distributed load (SN003a Table 3.2)', route:'uniform', c1in};
      if(shape.central && Math.abs(r-0.50)<=0.02) return {C1:1.348, C2:0.630, label:'simply supported + central point load (SN003a Table 3.2)', route:'point', c1in};
    }
    if(endFixed){
      if(shape.uniform) return {C1:2.578, C2:1.554, label:'fixed-ended + uniformly distributed load (SN003a Table 3.2)', route:'fixed-uniform', c1in};
      if(shape.central) return {C1:1.683, C2:1.645, label:'fixed-ended + central point load (SN003a Table 3.2)', route:'fixed-point', c1in};
    }
  }
  { // general diagram: Serna et al. quarter-point expression (SCI, NSC Nov 2013)
    const M2=Mst(xa+Ls/4), M3=Mst(xa+Ls/2), M4=Mst(xa+3*Ls/4);
    const c=sernaC1(Math.max(Mm,1e-9),M2,M3,M4);
    return {C1:c, C2:null, label:'general moment diagram &mdash; Serna et al. quarter-point expression (SCI): M(L/4)='+M2.toFixed(1)+', M(L/2)='+M3.toFixed(1)+', M(3L/4)='+M4.toFixed(1)+', M<sub>max</sub>='+Mm.toFixed(1)+' kN&middot;m', route:'serna', c1in};
  }
}
// ---- NCCI SN006a-EN-EU cantilever Mcr (doubly symmetric I/H), pure ----
// Mcr = C * Mcr0, Mcr0 = (pi/L) sqrt(E Iz G It); C from Tables 3.1-3.3 via
// sn006C() with kwt and eta; q + F combined by Eq (7). The SN006a boundary
// conditions replace the effective-length machinery (LE factor and the
// destabilising switch are NOT applied; load height enters through eta,
// warping through the root condition). `factors` = the combination's load
// factors (default a.governM), `zg` = the load height in mm for eta (default
// stdZgFor(factors)). Returns Mcr = null with `reason` (HTML) when the loading
// or the table range is not covered.
function mcrSN006aFor(a,sec,factors,zg){
  const E=a.E, G=81000, Iz=sec.Iy*1e4, It=sec.J*1e4, Iw=(sec.Iw||0)*1e12;
  const Lc=a.L;
  const Mcr0=Math.PI/Lc*Math.sqrt(E*Iz*G*It);        // N.mm
  const kwt=Math.sqrt(E*Iw/(G*It))/Lc;
  const hs=sec.D-sec.tf;
  const gcb=asCombo(factors||a.governM.combo), gfac=gcb.factors;
  const za=(zg!=null)? +zg : stdZgFor(gcb).zg;
  const eta=za/(hs/2);
  const warp=(endsList()[0].warp)?'restr':'free';   // root warping condition = End 1 warping flag
  // classify tip loading from the loads (2% de-minimis on the support moment)
  let Mq=0,MF=0,Mm2=0,nF2=0,nM2=0,other=false;
  const Lmm=S.L*1000;
  comboLoadPieces(gcb).forEach(p=>{
    const f=p.factor; if(!f) return; // zero-factor loads do not shape this combination
    if(p.type==='udl'&&p.x1<=1e-3&&Math.abs(p.x2-Lmm)<=1e-3) Mq+=p.w1*f*S.L*S.L/2;
    else if(p.type==='point'&&Math.abs(p.pos-Lmm)<=0.02*Lmm){ MF+=p.P*f*S.L; nF2++; }
    else if(p.type==='moment'&&Math.abs(p.pos-Lmm)<=0.02*Lmm){ Mm2+=Math.abs(p.M)*f; nM2++; }
    else other=true;
  });
  Mq+=(a.swPerM||0)*(gfac.G??0)*S.L*S.L/2;
  const tot=Math.abs(Mq)+Math.abs(MF)+Mm2;
  const dm=0.02*Math.max(tot,1e-9);
  const hasQ=Math.abs(Mq)>dm, hasF=Math.abs(MF)>dm, hasM=Mm2>dm;
  let C=null,Cq=null,CF=null,caseLbl='';
  if(other||hasM&&(hasQ||hasF)||nM2>1){ C=null; caseLbl='loading outside SN006a Tables 3.1-3.3'; }
  else if(hasM&&!hasQ&&!hasF){ C=sn006C('M',warp,kwt,0); caseLbl='external moment at the free end (Table 3.3)'; }
  else if(hasQ&&hasF){ Cq=sn006C('q',warp,kwt,eta); CF=sn006C('F',warp,kwt,eta);
    if(Cq!=null&&CF!=null) C=(Math.abs(Mq)+Math.abs(MF))/(Math.abs(Mq)/Cq+Math.abs(MF)/CF);
    caseLbl='uniform load + tip point load, interaction Eq (7)'; }
  else if(hasQ){ C=sn006C('q',warp,kwt,eta); caseLbl='uniformly distributed load (Table 3.1)'; }
  else if(hasF){ C=sn006C('F',warp,kwt,eta); caseLbl='point load at the free end (Table 3.2)'; }
  else { C=sn006C('q',warp,kwt,eta); caseLbl='self-weight only (Table 3.1)'; }
  let reason=null;
  if(C==null){
    if(kwt>1) reason="Cantilever LTB: &kappa;<sub>wt</sub> = "+kwt.toFixed(2)+" exceeds the SN006a table range (0&ndash;1); use a longer cantilever, a torsionally stiffer section, or a specialist tool (specialist software).";
    else if(Math.abs(eta)>0&&(eta<-2||eta>3)) reason="Cantilever LTB: load-height parameter &eta; = "+eta.toFixed(2)+" is outside the SN006a table range (&minus;2 to +3).";
    else reason="Cantilever LTB: "+caseLbl+" &mdash; not covered by SN006a; PASS is blocked (conservative C1=1.0 route removed in favour of the published method).";
  }
  return {Mcr0,kwt,hs,eta,za,warp,caseLbl,C,Cq,CF,Mq,MF,Mcr:(C==null? null : C*Mcr0),reason};
}
// ---- End conditions the closed-form Mcr chain can describe (pure, reads S.ends) ----
// The SN003a form (and the channel kappa chain and the box chain with I_w = 0)
// takes k = k_w = 1: FORK ends - U_y (lateral translation) and R_x (twist)
// held - at BOTH ends. NCCI SN006a describes the cantilever of
// isSn006aCantilever() (root U_y + R_z + R_x, free tip), doubly symmetric
// I/H only. A square hollow section is not susceptible to LTB whatever its
// ends (cl 6.3.2.1(2)). Every other end flag set is refused on the standard
// route (19 Sep 2026 review finding F-C: an end releasing R_x or U_y was
// given the fork-fork closed form, 0.41-0.68 of the eigenvalue), the printed
// chain still assuming fork ends; the FE eigen route models the flags.
// Laterally clamped (R_z) and warping-fixed ends raise M_cr; taking them as
// forks is conservative and is stated (the SN003a k = 0.5 rows are not
// applied). Returns {ok, sn006, msg (HTML, blocking when !ok), note (HTML,
// advisory), ends (HTML list of the released DOFs)}.
function stdMcrEndsStatus(st,sec){
  st=st||S;
  const ends=endsList(st);
  const sym={uy:'U<sub>y</sub>',rx:'R<sub>x</sub>',rz:'R<sub>z</sub>',warp:'warping'};
  const releases=e=>['uy','rx'].filter(k=>!e[k]).map(k=>sym[k]);
  const holds=e=>['uy','rx','rz','warp'].filter(k=>e[k]).map(k=>sym[k]);
  const cant=isCantilever(st);
  const useMcr='Use the FE eigenvalue M<sub>cr</sub> method (Axial &amp; lateral-torsional buckling), which models the end degrees of freedom; PASS is blocked on the standard route.';
  if(cant && sec && sec.kind==='I'){
    if(isSn006aCantilever(st)) return {ok:true, sn006:true, msg:null, note:null, ends:''};
    const [e1,e2]=ends;
    const rootTxt=releases(e1).length? 'End 1 (root) releases '+releases(e1).join(', ')+(e1.rz? '' : (releases(e1).length? ' and ' : 'releases ')+sym.rz) : (e1.rz? '' : 'End 1 (root) releases '+sym.rz);
    const tipTxt=holds(e2).length? 'End 2 (tip) restrains '+holds(e2).join(', ') : '';
    const desc=[rootTxt,tipTxt].filter(Boolean).join('; ');
    return {ok:false, sn006:false, kind:'sn006', ends:desc,
      msg:'Cantilever LTB on the standard route: NCCI SN006a (Tables 3.1&ndash;3.3) describes a cantilever whose root holds U<sub>y</sub>, R<sub>z</sub> and R<sub>x</sub> (v = v&prime; = &phi; = 0, warping restrained or free) and whose tip is free of every lateral restraint; here '+desc+', which the tables do not cover. '+useMcr,
      note:null};
  }
  if(sec && sec.isBox && st.family==='shs') return {ok:true, sn006:false, msg:null, note:null, ends:''};
  const bad=ends.filter(e=>!(e.uy&&e.rx));
  if(bad.length){
    const desc=bad.map(e=>'End '+e.n+' releases '+releases(e).join(' and ')).join('; ');
    return {ok:false, sn006:false, kind:'fork', ends:desc,
      msg:'Standard (closed-form) M<sub>cr</sub> for LTB: the SN003a form'+(sec&&sec.kind==='channel'? ' and the P385/P362 channel &kappa; chain assume' : ' assumes')+' fork ends &mdash; U<sub>y</sub> (lateral translation) and R<sub>x</sub> (twist) held &mdash; at both ends (k = k<sub>w</sub> = 1); here '+desc+'. The chain printed below assumes fork ends and would be unconservative for these end conditions. '+useMcr,
      note:null};
  }
  // a warping flag on an I_w = 0 box is not a boundary condition of the twist equation (warpingApplies, js/03-state-ui.js)
  const warpOn = !sec || warpingApplies(sec);
  const clamped=ends.filter(e=>e.rz).map(e=>'End '+e.n), warped=ends.filter(e=>e.warp&&warpOn).map(e=>'End '+e.n);
  let note=null;
  if(clamped.length||warped.length){
    const parts=[];
    if(clamped.length) parts.push('laterally clamped end'+(clamped.length>1?'s':'')+' (R<sub>z</sub> held at '+clamped.join(' and ')+')');
    if(warped.length) parts.push('warping-fixed end'+(warped.length>1?'s':'')+' ('+warped.join(' and ')+')');
    note='Standard (closed-form) M<sub>cr</sub>: '+parts.join(' and ')+' taken as fork end'+(clamped.length+warped.length>1?'s':'')+', k = k<sub>w</sub> = 1 (conservative: the SN003a k = 0.5 / k<sub>w</sub> = 0.5 rows are not applied; the FE eigen route uses v&prime; = 0 / &phi;&prime; = 0 and gives the higher M<sub>cr</sub>).';
  }
  const warpIgnored=ends.filter(e=>e.warp&&!warpOn).map(e=>'End '+e.n);
  if(warpIgnored.length) note=(note? note+' ' : '')+'Warping flag at '+warpIgnored.join(' and ')+' not applied: this closed section has I<sub>w</sub> = 0, so &phi;&prime; is not a boundary condition of its twist equation (EN 1993-1-1 6.2.7(7)).';
  return {ok:true, sn006:false, msg:null, note, ends:''};
}
// ---- Standard (closed-form) Mcr for one segment, pure ----
// Used by the eigen method for the "Mcr eigen / Mcr standard" comparison (no
// eigen solve needed) and by the report. Cantilever (I/H) -> SN006a; otherwise
// the SN003a form with C1 from sn003aC1 over the segment, LE = LE-factor
// (x1.2 if destabilising) x segment length, z_g = the load height of the
// combination (stdZgFor) with C2 where published. End flag sets the closed
// form cannot describe (stdMcrEndsStatus) return Mcr = null with the reason.
// `diag` = {fb, factors} selects the combination whose diagram and loads the
// value describes (default a.governM), so the eigen method can pass its own
// LTB-governing combination and get a sign-consistent comparison.
// Returns kN.m: {route, Mcr (null if not covered), C1, C2, label, c1in, seg,
// LE, zg, zgSource, zgUsed, zgNote (HTML), zgBlocked}.
function mcrStandardFor(a,sec,seg,diag){
  const isCant=isCantilever(S);
  seg=seg||c1Segment(a);
  const fb=(diag&&diag.fb)||a.governM.fb;
  const fac=(diag&&diag.combo)||((diag&&diag.factors)? {factors:diag.factors} : a.governM.combo);   // combination (pattern-aware)
  const zgi=stdZgFor(fac);
  const es=stdMcrEndsStatus(S,sec);
  if(!es.ok){
    return {route:'unsupported', Mcr:null, C1:null, C2:null, label:'closed form not applicable to these end conditions ('+es.ends+'; '+(es.kind==='sn006'? 'SN006a needs a root holding U<sub>y</sub>, R<sub>z</sub>, R<sub>x</sub> and a free tip' : 'the SN003a form needs fork ends U<sub>y</sub> + R<sub>x</sub> at both ends')+')',
      c1in:c1Inputs(fb,seg.xa,seg.xb), seg, LE:null, endsMsg:es.msg, zg:zgi.zg, zgSource:zgi.source, zgUsed:false, zgNote:'', zgBlocked:false};
  }
  if(es.sn006){
    const r=mcrSN006aFor(a,sec,fac,zgi.zg);
    const c1in=c1Inputs(fb,seg.xa,seg.xb);
    return {route:'sn006a', Mcr:(r.Mcr!=null? r.Mcr/1e6 : null), C1:r.C, C2:null,
      label:'cantilever SN006a &mdash; '+r.caseLbl, c1in, seg, LE:a.L, sn006:r,
      zg:zgi.zg, zgSource:zgi.source, zgUsed:Math.abs(zgi.zg)>1e-9, zgNote:'z<sub>g</sub> = '+(zgi.zg>0?'+':'')+zgi.zg.toFixed(0)+' mm through &eta; = '+r.eta.toFixed(2)+' (SN006a)', zgBlocked:false};
  }
  const c1r=sn003aC1(a,isCant,seg.whole? undefined : seg,{fb,combo:fac});
  const LE=ltbLeFactor()*(S.destab?1.2:1)*(seg.xb-seg.xa);
  const cf=mcrClosedForm(sec,a.E,LE,c1r.C1,c1r.C2,zgi.zg);
  const zs=stdZgStatus(zgi.zg,c1r.C2,cf.zgUsed);
  const route= sec.kind==='channel'? 'channel' : c1r.route;
  const endsTag= es.note? ' [k = k<sub>w</sub> = 1: clamped / warping-fixed end(s) taken as forks]' : '';
  return {route, Mcr:cf.Mcr/1e6, C1:c1r.C1, C2:c1r.C2, label:c1r.label+endsTag, c1in:c1r.c1in, seg, LE,
    T1:cf.T1/1e3, zg:zgi.zg, zgSource:zgi.source, zgUsed:cf.zgUsed, zgNote:zs.text, zgBlocked:!!zs.block, endsNote:es.note};
}

function checksEC3UnrestrainedSCI(a){
  // STANDARD METHOD (closed-form Mcr). Kept alive by js/08-mcr-eigen-patch.js
  // as window.checksEC3UnrestrainedStandard and used when S.mcrMethod ===
  // 'standard'; the patch's own function (FE eigensolver) is the default.
  // SCI worked-example procedure: unrestrained beam to BS EN 1993-1-1 (UK NA).
  // Cross-section checks are identical to the restrained case. LTB design
  // basis for I/H sections: the elastic critical moment route, Mcr from the
  // SN003a closed form (k = kw = 1, G = 81000 N/mm2, C2*zg where published),
  // lamLT = sqrt(Wy fy / Mcr), chiLT per cl 6.3.2.3 with lamLT0 = 0.4,
  // beta = 0.75 (NA 2.17), curve from NA Table NA.1 (ltbCurveNA), then
  // chiLT,mod = chiLT/f with kc = 1/sqrt(C1) (NA 2.18) and Mb,Rd <= Mc,Rd.
  // The P362 Expn (6.55) simplified slenderness lamLT = (1/sqrt(C1)) 0.9
  // lamZbar sqrt(betaW) is evaluated alongside as a COMPARISON only (ltb.MbSimp);
  // it never sets the verdict. Boxes take the same Mcr chain with Iw = 0,
  // cantilevers SN006a, channels the P385/P362 kappa chain.
  const b=checksEC3Restrained(a);
  const sec=a.sec, fy=a.py, E=a.E, gM1=1.0;
  const unsupported=b.unsupported.slice();
  const advisory=(b.advisory||[]).slice();
  const Wy=b.Wy;
  const isCant=isCantilever(S);
  // End conditions the closed form can describe (19 Sep 2026 review finding
  // F-C): fork ends U_y + R_x at both ends, or the SN006a cantilever of
  // isSn006aCantilever() for an I/H section; anything else is NOT VERIFIED
  // on this route (message below, the fork-ended chain still printed), a
  // clamped / warping-fixed end is taken as a fork and said so.
  const endsStd=stdMcrEndsStatus(S,sec);
  const LE=ltbLeFactor()*(S.destab?1.2:1)*a.L;
  const gfac=a.governM.combo;     // governing-moment combination (pattern-aware load list)
  const c1r=sn003aC1(a,isCant);   // whole member: the closed form treats the member as one segment
  const C1=c1r.C1, kcr=kcFromC1(C1), invSqrtC1=kcr.kcRaw, kc=kcr.kc;   // k_c floored at 1/sqrt(2.76) = 0.60 (Table 6.6 lower bound)
  const zgi=stdZgFor(gfac);       // load height of the governing combination (per-load z_g, most destabilising)
  const zgStd=zgi.zg;
  const segStd=c1Segment(a);
  if(!segStd.whole) advisory.push("Standard (closed-form) M<sub>cr</sub>: the member is treated as ONE segment of length L<sub>E</sub> = "+(ltbLeFactor()*(S.destab?1.2:1)).toFixed(2)+"&times;L with C<sub>1</sub> from the whole-member moment diagram; intermediate lateral restraints are not applied on this route (conservative on L<sub>E</sub>). Use the FE eigenvalue method for the bay-by-bay M<sub>cr</sub>, or verify each bay separately.");
  const curve=ltbCurveNA(sec);
  // chi_LT chain of cl 6.3.2.3 (NA 2.17: lamLT0 = 0.4, beta = 0.75) with the
  // NA 2.18 f-factor from kc = 1/sqrt(C1); f = 1 when `noF` (cantilevers)
  const chiChain=(lamLT,noF)=>{
    if(lamLT<=0.4) return {Phi:null,chi:1,f:1,chiMod:1,ign:true};
    const Phi=0.5*(1+curve.alphaLT*(lamLT-0.4)+0.75*lamLT*lamLT);
    let chi=1/(Phi+Math.sqrt(Math.max(Phi*Phi-0.75*lamLT*lamLT,1e-12)));
    chi=Math.min(chi,1,1/(lamLT*lamLT));
    let f=noF? 1 : Math.min(1-0.5*(1-kc)*(1-2*Math.pow(lamLT-0.8,2)),1);
    const chiMod=Math.min(chi/f,1,1/(lamLT*lamLT));
    return {Phi,chi,f,chiMod,ign:false};
  };
  let ltb, zgs=null;
  if(sec.isBox){
    // Closed section: the same SN003a chain with Iw = 0 (warping neglected),
    // C2*zg where published, chi_LT from the NA Table NA.1 curve (hot-finished
    // hollow sections share the I/H h/b allocation, cold-formed c/d - the same
    // ltbCurveNA() the eigen method uses). lamLT <= 0.4 reproduces the
    // exemption (cl 6.3.2.2(4); square hollow sections are also exempt by
    // cl 6.3.2.1(2)); a slender RHS on a long span gets a real chi_LT.
    const cf=mcrClosedForm(sec,E,LE,C1,c1r.C2,zgStd);
    const Mcr=cf.Mcr;
    const lamLT=Math.sqrt(Wy*fy/Mcr);
    const sB=chiChain(lamLT,false);
    const Mb=Math.min(sB.chiMod*Wy*fy/gM1/1e6,b.McRd);
    zgs=stdZgStatus(zgStd,c1r.C2,cf.zgUsed);
    ltb={na:sB.ign, box:true, T1:cf.T1/1e3, GIt:cf.GIt/1e9, Mcr:Mcr/1e6, lamLTmcr:lamLT, ignM:sB.ign,
      PhiM:sB.Phi, chiM:sB.chi, fM:sB.f, chiModM:sB.chiMod, curve, kc, kcRaw:kcr.kcRaw, kcFloored:kcr.floored, invSqrtC1,
      zg:zgStd, C2:c1r.C2, zgUsed:cf.zgUsed, MbSimp:Mb, MbMcr:Mb, MbRd:Mb};
  } else if(endsStd.sn006){
    // NCCI SN006a-EN-EU cantilever path (doubly symmetric I/H, root U_y + R_z +
    // R_x, free tip - isSn006aCantilever): see mcrSN006aFor(). An I/H
    // cantilever with other lateral flags falls through to the fork-ended
    // chain below and is refused by endsStd (NOT VERIFIED).
    const sn=mcrSN006aFor(a,sec,gfac,zgStd);
    const {Mcr0,kwt,eta,warp,caseLbl,C,Cq,CF,Mq,MF}=sn;
    zgs={text:'z<sub>g</sub> = '+(zgStd>0?'+':'')+zgStd.toFixed(0)+' mm through &eta; = z<sub>g</sub>/(h<sub>s</sub>/2) = '+eta.toFixed(2)+' (SN006a)', block:null, advisory:null};
    if(C==null){
      unsupported.push(sn.reason);
      ltb={na:false,cant:true,Mcr0:Mcr0/1e6,kwt,eta,warp,caseLbl,C:0,Cq,CF,Mq,MF,
        lamLTsimp:0,lamLTmcr:0,ignS:true,ignM:true,chiM:1,fM:1,chiModM:1,PhiM:null,
        curve,kc:1,invSqrtC1:1,MbSimp:0,MbMcr:0,MbRd:0,Mcr:0,T1:0,IwIz:0,GIt:0,zg:zgStd,zgUsed:Math.abs(zgStd)>1e-9};
    } else {
      const Mcr=C*Mcr0;
      const lamLT=Math.sqrt(Wy*fy/Mcr);
      const sB=chiChain(lamLT,true);   // kc/f-factor not applied for cantilevers (no published kc)
      const Mb=Math.min(sB.chiMod*Wy*fy/gM1/1e6,b.McRd);
      ltb={na:false,cant:true,Mcr0:Mcr0/1e6,kwt,eta,warp,caseLbl,C,Cq,CF,Mq,MF,
        lamLTsimp:lamLT,lamLTmcr:lamLT,PhiS:sB.Phi,chiS:sB.chi,ignS:sB.ign,
        PhiM:sB.Phi,chiM:sB.chi,fM:1,chiModM:sB.chiMod,ignM:sB.ign,curve,kc:1,invSqrtC1:1,
        lamZ:0,lam1:0,lamZbar:0,rootBw:1,hb:sec.D/sec.B,
        MbSimp:Mb,MbMcr:Mb,MbRd:Mb,Mcr:Mcr/1e6,McrBack:Mcr/1e6,T1:0,IwIz:0,GIt:0,fS:1,chiModS:sB.chiMod,zg:zgStd,zgUsed:Math.abs(zgStd)>1e-9};
    }
  } else if(sec.kind==='channel'){
    // P385 / P362 channel chain (spec 13.4): lamLT = (L/i_z)/kappa, curve d,
    // rolled-section chi formulas (validated: Ex 4 chi = 0.29), no f-factor;
    // M_cr back-calculated = W_y*f_y / lamLT^2 (needed for the Annex A k_alpha).
    const ry=sec.ry*10;
    const kappa=({S275:96,S355:85,S460:74})[S.grade]||96;
    const lamLT=(LE/ry)/kappa;
    const alphaLT=0.76; // curve d (non-doubly-symmetric)
    let chi=1,Phi=null,ign=true;
    if(lamLT>0.4){
      Phi=0.5*(1+alphaLT*(lamLT-0.4)+0.75*lamLT*lamLT);
      chi=Math.min(1/(Phi+Math.sqrt(Math.max(Phi*Phi-0.75*lamLT*lamLT,1e-12))),1,1/(lamLT*lamLT));
      ign=false;
    }
    const Mb=Math.min(chi*Wy*fy/gM1/1e6,b.McRd);
    const McrBack=Wy*fy/(lamLT*lamLT)/1e6;
    // Mcr route for channels: a PFC bent about its major axis is symmetric about
    // the axis of bending (Wagner term zj = 0), so with the load through the
    // SHEAR CENTRE and fork supports at both ends the doubly-symmetric Mcr
    // expression is theoretically exact (validated against an independent commercial-software
    // channel example: Mcr 330.9, chi 0.881, f 0.838, Mb = Mc 97.625). Route is
    // OFFERED only under those conditions: no active torsion (e = 0 everywhere),
    // z_g = 0, two end fork supports (not a cantilever). The kappa chain stays
    // the primary basis; the Mcr route rescues it, mirroring the old I-section pattern.
    let chanMcr=null;
    const chanMcrOK = !(a.tors&&a.tors.on) && Math.abs(zgStd)<1e-9 && !isCant &&
      endsList().every(e=>e.uz&&e.uy&&e.rx) && (sec.Iw||0)>0;   // vertically supported fork ends at both ends
    let MbMcr2=Mb;
    if(chanMcrOK){
      const G=81000, Iz=sec.Iy*1e4, It=sec.J*1e4, Iw=(sec.Iw||0)*1e12;
      const T1c=Math.PI*Math.PI*E*Iz/(LE*LE);
      const McrC=C1*T1c*Math.sqrt(Math.max(Iw/Iz+G*It/T1c,0));
      const lamC=Math.sqrt(Wy*fy/McrC);
      let PhiC=null,chiC=1,fC=1,chiModC=1,ignC=true;
      if(lamC>0.4){
        PhiC=0.5*(1+0.76*(lamC-0.4)+0.75*lamC*lamC);
        chiC=Math.min(1/(PhiC+Math.sqrt(Math.max(PhiC*PhiC-0.75*lamC*lamC,1e-12))),1,1/(lamC*lamC));
        fC=Math.min(1-0.5*(1-kc)*(1-2*Math.pow(lamC-0.8,2)),1);
        chiModC=Math.min(chiC/fC,1,1/(lamC*lamC)); ignC=false;
      }
      MbMcr2=Math.min(chiModC*Wy*fy/gM1/1e6, b.McRd);
      chanMcr={Mcr:McrC/1e6,lam:lamC,Phi:PhiC,chi:chiC,f:fC,chiMod:chiModC,ign:ignC,Mb:MbMcr2,T1:T1c/1e3};
    }
    zgs={text: Math.abs(zgStd)<1e-9? 'z<sub>g</sub> = 0 (load through the shear centre)' : 'z<sub>g</sub> = '+(zgStd>0?'+':'')+zgStd.toFixed(0)+' mm entered &mdash; not used by the P385/P362 &kappa; chain (load height enters the channel route only through the destabilising L<sub>E</sub> switch); the shear-centre M<sub>cr</sub> route is not offered', block:null, advisory:null};
    ltb={na:false,channel:true,chanMcr,ry,kappa,lamLTsimp:lamLT,lamLTmcr:lamLT,PhiS:Phi,chiS:chi,ignS:ign,
      PhiM:Phi,chiM:chi,fM:1,chiModM:chi,ignM:ign,curve:{alphaLT:0.76,curve:'d'},kc,kcRaw:kcr.kcRaw,kcFloored:kcr.floored,invSqrtC1,
      lamZ:LE/ry,lam1:0,lamZbar:0,rootBw:1,hb:sec.D/sec.B,
      MbSimp:Mb,MbMcr:MbMcr2,MbRd:Mb,Mcr:McrBack,McrBack,T1:0,IwIz:0,GIt:0,fS:1,chiModS:chi,zg:zgStd,zgUsed:false};
  } else {
    const ry=sec.ry*10;
    const lamZ=LE/ry;
    const lam1=Math.PI*Math.sqrt(E/fy);
    const lamZbar=lamZ/lam1;
    const rootBw=b.cl.cls<=2? 1 : Math.sqrt(sec.Zx/sec.Sx);
    const hb=sec.D/sec.B;
    // (B) DESIGN BASIS: elastic critical moment (SN003a; k = kw = 1; C2*zg
    // load-height term where published for the recognised diagram) - mcrClosedForm()
    const C2=c1r.C2;
    const cf=mcrClosedForm(sec,E,LE,C1,C2,zgStd);
    const T1=cf.T1, IwIz=cf.IwIz, GIt=cf.GIt, zgUsed=cf.zgUsed, Mcr=cf.Mcr; // N, mm2, N.mm2, N.mm
    const lamLTmcr=Math.sqrt(Wy*fy/Mcr);
    const sB=chiChain(lamLTmcr,false);
    const MbMcr=Math.min(sB.chiMod*Wy*fy/gM1/1e6,b.McRd);
    // (A) COMPARISON ONLY: simplified slenderness, P362 Expn (6.55)
    const lamLTsimp=invSqrtC1*0.9*lamZbar*rootBw;
    const sA=chiChain(lamLTsimp,false);
    const MbSimp=Math.min(sA.chiMod*Wy*fy/gM1/1e6,b.McRd);
    zgs=stdZgStatus(zgStd,C2,zgUsed);
    ltb={na:false,ry,lamZ,lam1,lamZbar,rootBw,hb,curve,kc,kcRaw:kcr.kcRaw,kcFloored:kcr.floored,invSqrtC1,
      lamLTsimp,PhiS:sA.Phi,chiS:sA.chi,fS:sA.f,chiModS:sA.chiMod,ignS:sA.ign,MbSimp,
      T1:T1/1e3,IwIz:IwIz/1e2,GIt:GIt/1e9,Mcr:Mcr/1e6,zg:zgStd,C2,zgUsed,
      lamLTmcr,PhiM:sB.Phi,chiM:sB.chi,fM:sB.f,chiModM:sB.chiMod,ignM:sB.ign,MbMcr,
      MbRd:MbMcr};
  }
  // ---- load height bookkeeping (always printed; may block PASS) ----
  ltb.zgSource=zgi.source;
  ltb.zgNote=zgs? zgs.text : '';
  ltb.zgBlocked=!!(zgs&&zgs.block);
  if(zgs&&zgs.block) unsupported.push(zgs.block);
  if(zgs&&zgs.advisory) advisory.push(zgs.advisory);
  // ---- end conditions of the closed form (always recorded; may block PASS) ----
  ltb.endsOk=endsStd.ok; ltb.endsMsg=endsStd.msg; ltb.endsNote=endsStd.note;
  if(!endsStd.ok) unsupported.push(endsStd.msg);
  if(endsStd.note) advisory.push(endsStd.note);
  // ---- method tags and the reference software-style C1 inputs (M1, M2, Mo, psi, mu) ----
  // The closed form derives C1 from the WHOLE member, so the printed inputs are
  // the whole-member values (segment 0..L); the eigen method fills the same
  // fields for its governing span. McrEigen is null here: the eigensolver is
  // not run on the standard route (keeps it fast).
  ltb.mcrMethod='standard';
  ltb.c1in=c1r.c1in;
  ltb.c1seg={xa:0,xb:a.L,whole:true};
  ltb.c1route= ltb.cant? 'sn006a' : ltb.channel? 'channel' : c1r.route;
  ltb.c1label= ltb.cant? ('cantilever SN006a &mdash; '+ltb.caseLbl) : ltb.channel? ('channel &mdash; P385/P362 kappa chain; '+c1r.label) : c1r.label;
  ltb.C1show= ltb.cant? (ltb.C||0) : C1;
  ltb.McrStandard= ltb.cant? (ltb.Mcr>0? ltb.Mcr : null) : ltb.channel? (ltb.chanMcr? ltb.chanMcr.Mcr : ltb.McrBack) : ltb.Mcr;
  ltb.McrEigen=null;
  ltb.segWhole=segStd.whole;
  const Mx=b.Mx;
  let ltbUtil = ltb.MbRd>0? Mx/ltb.MbRd : 0;
  let ltbBasis;
  if(ltb.box){
    ltbBasis = ltb.ignM
      ? 'closed hollow section, &lambda;&#772;<sub>LT</sub> = '+ltb.lamLTmcr.toFixed(3)+' &le; 0.4 from the SN003a M<sub>cr</sub> with I<sub>w</sub> = 0: lateral-torsional buckling effects may be ignored (cl 6.3.2.2(4)'+(S.family==='shs'? '; a square hollow section is also not susceptible by cl 6.3.2.1(2)' : '')+'), M<sub>b,Rd</sub> = M<sub>c,Rd</sub>'
      : 'closed hollow section on a long span: M<sub>cr</sub> method (SN003a with I<sub>w</sub> = 0), &lambda;&#772;<sub>LT</sub> = '+ltb.lamLTmcr.toFixed(3)+' &gt; 0.4, &chi;<sub>LT,mod</sub> from cl 6.3.2.3 with curve '+curve.curve+' (NA Table NA.1, hot-finished / cold-formed hollow section) and k<sub>c</sub> = 1/&radic;C<sub>1</sub>; M<sub>b,Rd</sub> = &chi;<sub>LT,mod</sub>W<sub>y</sub>f<sub>y</sub>/&gamma;<sub>M1</sub> &le; M<sub>c,Rd</sub>';
  } else if(ltb.cant){
    ltbBasis = 'M<sub>cr</sub> method (NCCI SN006a cantilever): M<sub>cr</sub> = C&middot;M<sub>cr,0</sub>, &lambda;&#772;<sub>LT</sub> = &radic;(W<sub>y</sub>f<sub>y</sub>/M<sub>cr</sub>), &chi;<sub>LT</sub> from cl 6.3.2.3 curve '+curve.curve+' with f = 1; M<sub>b,Rd</sub> = &chi;<sub>LT</sub>W<sub>y</sub>f<sub>y</sub>/&gamma;<sub>M1</sub> &le; M<sub>c,Rd</sub>';
  } else if(ltb.channel){
    ltbBasis = 'P385/P362 channel &kappa; chain, &lambda;&#772;<sub>LT</sub> = (L<sub>E</sub>/i<sub>z</sub>)/&kappa;, curve d, no f-factor; M<sub>b,Rd</sub> = &chi;<sub>LT</sub>W<sub>y</sub>f<sub>y</sub>/&gamma;<sub>M1</sub> &le; M<sub>c,Rd</sub>';
    if(ltbUtil>1.0001 && ltb.MbMcr>0 && Mx/ltb.MbMcr<=1.0001){
      ltbUtil=Mx/ltb.MbMcr; ltb.MbRd=ltb.MbMcr;
      ltbBasis='M<sub>cr</sub> method (SN003a with z<sub>j</sub> = 0, load through the shear centre) &mdash; the P385/P362 channel &kappa; chain is exceeded, but it is conservative; adequacy is demonstrated by the M<sub>cr</sub> route';
    }
  } else {
    ltbBasis = 'M<sub>cr</sub> method (SN003a closed form, k = k<sub>w</sub> = 1, G = 81000 N/mm&sup2;'+(ltb.zgUsed? ', C<sub>2</sub>z<sub>g</sub> load-height term' : '')+'): &lambda;&#772;<sub>LT</sub> = &radic;(W<sub>y</sub>f<sub>y</sub>/M<sub>cr</sub>) = '+ltb.lamLTmcr.toFixed(3)+', &chi;<sub>LT,mod</sub> from cl 6.3.2.3 curve '+curve.curve+' with k<sub>c</sub> = 1/&radic;C<sub>1</sub> (NA 2.18); M<sub>b,Rd</sub> = &chi;<sub>LT,mod</sub>W<sub>y</sub>f<sub>y</sub>/&gamma;<sub>M1</sub> &le; M<sub>c,Rd</sub>. The P362 Expn 6.55 simplified slenderness (&lambda;&#772;<sub>LT</sub> = '+ltb.lamLTsimp.toFixed(3)+', M<sub>b,Rd</sub> = '+ltb.MbSimp.toFixed(1)+' kN&middot;m) is printed for comparison only';
  }
  // end conditions of the closed form: fork ends assumed (a released end is NOT VERIFIED), clamped / warping-fixed ends taken as forks
  if(!endsStd.ok) ltbBasis+='. NOT VERIFIED: the closed form assumes fork ends (U<sub>y</sub> + R<sub>x</sub> at both ends'+(isCant? ', or the SN006a cantilever root / free tip' : '')+') and '+endsStd.ends+' &mdash; use the FE eigenvalue method';
  else if(endsStd.note) ltbBasis+='. Clamped / warping-fixed end(s) taken as fork ends, k = k<sub>w</sub> = 1 (conservative)';
  // ---- BS EN 1993-6 Annex A: LTB + minor-axis bending + torsion interaction (P385 6.2/8.2) ----
  let annex=null;
  if(b.tor && b.tor.p385 && !ltb.na){
    // chi_LT WITHOUT the f-factor (P385 validation basis), curve per table
    const chiA=ltb.channel? ltb.chiM : (ltb.ignM? 1 : ltb.chiM);
    const MbA=chiA*Wy*fy/gM1/1e6;
    const McrA=ltb.channel? ltb.McrBack : ltb.Mcr;
    // Cmz proxy for the minor-axis diagram: the tabulated load case behind C1
    // (central point load 0.9, UDL 0.95, otherwise 1.0 conservative)
    const CmzProxy = (c1r.route==='point')? 0.9 : (c1r.route==='uniform')? 0.95 : 1.0;
    // resistances CLASS-CONSISTENT (elastic for Class 3) - required: with plastic
    // values a Class 3 member scores unconservatively (exposed by an independent commercial-software
    // UC 152x152x23 warping-torsion example: elastic gives 1.03-1.05 FAIL, plastic
    // would have shown 0.91 PASS).
    // 20 Sep 2026 torsion + N/Mz: Eq (A.1) per station with M_z,Ed = M_z,tot = imposed M_z + phi.M_y
    // in both the C_mz M_z/M_z,Rd term and k_zw; C_mz = 1.0 when an M_z is imposed (annexAEval above)
    annex=annexAEval(b.tor,MbA,McrA,CmzProxy);
    // 20 Sep 2026: M_y,Ed >= M_cr is a FAILURE (LTB governs, the Annex A utilisation is carried as 99), not an
    // unverified check - reported through advisory with the FAIL: prefix so the brief prints a Warning row
    if(annex.unbounded) advisory.push(KALPHA_UNBOUNDED_FAIL);
  }
  // member buckling: susceptible to torsional deformation unless closed section
  // or LTB plays no part (chiLT = 1); cantilever/channel handled per path.
  // The Mb,Rd handed to Eq 6.61/6.62 is the design value above (Mcr route).
  const isCantU=isCantilever(S);
  const useB1u = sec.isBox || ltb.na || (ltb.MbRd>=b.McRd*0.9999);
  const buck=(b.ax && !b.ax.tension)? annexB2(a,sec,fy,b.cl,ltb.MbRd>0? ltb.MbRd : b.McRd,useB1u,isCantU,b.aeff) : null;
  if(buck && buck.tfb && !buck.tfb.ok) unsupported.push('PFC under axial compression: '+buck.tfb.reason+'; torsional / torsional-flexural buckling (cl 6.3.1.4) cannot be verified, PASS is blocked.');
  const restraintF=restraintForces(a,sec);
  const utils=[
    {name:"Shear  V_Ed/V_c,Rd",val:b.shearUtil},
    {name:"Bending  M_Ed/M_c,Rd",val:b.momUtil},
    {name:"LTB  M_Ed/M_b,Rd",val:ltbUtil},
    {name:"Deflection",val:b.dmax/b.dlimit},
  ];
  if(b.ax){
    utils.push({name: b.ax.tension? "Tension  N_Ed/N_t,Rd" : (b.ax.aeff&&b.ax.aeff.active? "Compression  N_Ed/N_c,Rd (A_eff)" : "Compression  N_Ed/N_pl,Rd"), val:b.ax.nUtil});
    utils.push({name: b.ax.biax? ("Biaxial bending"+((S.axial||0)!==0?" + axial":"")+" (6.2.9.1)") : "Bending+axial cross-section (6.2.9)",val:b.ax.mUtil});
    // Eq 6.61/6.62 are needed with axial compression AND for biaxial bending on
    // an LTB-susceptible member with N_Ed = 0 (kzy*My/MbRd + kzz*Mz/MczRd).
    if(!b.ax.tension && buck && (buck.Fc>1e-9 || buck.biax)){
      if(buck.tfb && buck.tfb.ok) utils.push({name:TFB_UTIL_NAME,val:buck.tfb.util});
      utils.push({name:"Member buckling y-y (Eq 6.61)",val:buck.u1});
      utils.push({name:"Member buckling z-z (Eq 6.62)",val:buck.u2});
    }
  }
  if(annex) utils.push({name:"LTB+torsion (EN 1993-6 Annex A)",val:annex.u});
  if(b.coex) utils.push({name:b.coex.pureShearFail? "Pure shear failure at M-V check point (6.2.6)" : "Bending+shear coexistent (6.2.8)",val:b.coex.u});
  if(b.mvn) utils.push({name:mvnUtilName(b.mvn),val:b.mvn.u});
  if(b.web&&b.web.checked){
    utils.push({name:WEB_UTIL_NAMES[0],val:b.web.util2});
    utils.push({name:WEB_UTIL_NAMES[1],val:b.web.util72});
  }
  if(b.tor&&b.tor.box){
    utils.push({name:"Torsion  T_Ed/T_Rd",val:b.tor.torUtil});
    utils.push({name:"Shear+torsion  V_Ed/V_pl,T,Rd",val:b.tor.vtUtil});
  }
  if(b.tor&&b.tor.p385){
    utils.push({name:"Bending+torsion cross-section (P385 3.1.2)",val:b.tor.cross.u});
    utils.push({name:"Shear+torsion  V_Ed/V_pl,T,Rd",val:b.tor.vtUtil});
  }
  // 20 Sep 2026 torsion + N/Mz: elastic yield criterion (6.1), cl 6.2.7(5); 20 Sep 2026 review: in utils only when
  // verdict-binding (elasticBindingPolicy in checksEC3Restrained), otherwise in b.info (carried over) with its advisory
  if(b.tor&&b.tor.elastic&&b.tor.elastic.binding) utils.push({name:ELASTIC_TORSION_UTIL_NAME,val:b.tor.elastic.u});
  // 20 Sep 2026 torsion + N/Mz: advisory superposition of Eq 6.62 and (A.1) with this path's k_alpha (information only)
  if(b.tor&&b.tor.p385){ b.tor.superposition=torsionSuperposition(b.tor,annex,buck); if(b.tor.superposition) advisory.push(b.tor.superposition.text); }
  // 20 Sep 2026 review: the basis text of this path (6.61/6.62 and (A.1) as evaluated here)
  if(b.tor&&(b.tor.p385||b.tor.box)) b.tor.combinedBasis=torsionCombinedBasis({box:!!b.tor.box,tension:!!(b.ax&&b.ax.tension),hasN:Math.abs(S.axial||0)>1e-9,hasMz:Math.abs(S.Mz||0)>1e-9,restrained:false,buckEvaluated:!!(!(b.ax&&b.ax.tension)&&buck&&(buck.Fc>1e-9||buck.biax)),annexEvaluated:!!annex,binding:!!(b.tor.elastic&&b.tor.elastic.binding)});
  let gov=utils[0]; utils.forEach(u=>{ if(u.val>gov.val) gov=u; });
  const pass=unsupported.length===0 && utils.every(u=>u.val<=1.0001);
  return Object.assign({},b,{sci:false,sciU:true,mcrMethod:'standard',unsupported,advisory,ltb,ltbUtil,ltbBasis,C1,c1label:c1r.label,LE,utils,gov,pass,annex,buck,restraintForces:restraintF});
}


/* ==== beam-v03 module 12 ==== */
/* ===========================================================================
   Warping-torsion finite-element solver for open sections
   (19 Sep 2026 gap closure, group G4; docs/COVERAGE_MATRIX.md section 5 item
   11, trigger 2.6 of docs/EC3_BEAM_TRIGGER_LIST.md)
   ---------------------------------------------------------------------------
   Governing equation (Vlasov, linear elastic, small twist - the equation
   behind the SCI P385 Appendix C closed forms):

       E I_w phi''''(x) - G I_T phi''(x) = m_t(x)

   phi  = twist [rad], x along the member [mm], m_t = applied torque per unit
   length [N.mm/mm] (positive in the same sense as phi). Internal actions:
       T(x)   = G I_T phi' - E I_w phi'''   total torque   [N.mm]
       T_t(x) = G I_T phi'                  St Venant part [N.mm]
       B(x)   = E I_w phi''                 bimoment       [N.mm2]
       T_w(x) = -E I_w phi'''               warping part   [N.mm]
   so that dT/dx = -m_t (a distributed torque decreases T, a point torque P
   at x_p gives T(x_p+) = T(x_p-) - P).

   Discretisation: two-node elements, DOFs (phi, phi') per node, cubic Hermite
   shape functions - exactly the element that js/08-mcr-eigen-patch.js uses
   for the lateral-torsional eigenproblem. Its kBend(EI, Le) and kTors(GJ,
   Le) matrices are private to that file's closure, so they are duplicated
   here as the two small pure functions torsionFeKBend / torsionFeKTors
   (identical entries; hand-checked in tests/torsion-fe.test.cjs):
       K_e = kBend(E I_w, Le) + kTors(G I_T, Le)
   Consistent nodal torques of a linearly varying m_t over an element (same
   integrals as the bending solver's distributed-load vector).

   Boundary conditions (from the end degree-of-freedom flags, 19 Sep 2026
   single-span scope: an entry per end whose twist R_x is restrained):
       twist restrained      phi = 0 (fork / torsional restraint); one such
                             end = torsion cantilever
       warping free          phi'' = 0 (natural, nothing imposed) - default,
                             matching P385 fork ends
       warping fixed         phi' = 0 imposed (the end's warping flag)
       free end / tip        natural: T = 0 and B = 0 come out of the
                             assembled equilibrium exactly
       in-plane hinge        twist and warping continuous (a bending release
                             does not release twist)

   Stress recovery (per station):
       phi, phi'   nodal DOFs
       B           element end forces f = K_e d_e - f_eq,e (the consistent
                   flux; satisfies the natural conditions B = 0 at free and
                   fork ends exactly) averaged at interior nodes where B is
                   continuous. At an interior node whose phi' DOF is fixed
                   (warping restrained) the restraint applies a REACTION
                   BIMOMENT, so B jumps there exactly like T jumps at a point
                   torque: the node is reported as three stations with the
                   left element value at x - 0.01, the larger-magnitude side
                   at x and the right element value at x + 0.01 mm (averaging
                   across the jump under-reported the peak B_Ed and made the
                   mesh-convergence measure creep instead of converge;
                   19 Sep 2026 review finding)
       T           element end forces (nodal equilibrium exact); a node where
                   T jumps (point torque) is reported as
                   three stations x - 0.01, x, x + 0.01 mm like p385Solve
       phi''       B/(E I_w);  phi''' = (G I_T phi' - T)/(E I_w)
   The returned object has the shape of p385Solve's result ({xs, phi, p1, p2,
   p3, X}) so the P385 3.1.2 cross-section check, the V_pl,T,Rd sweep and
   the EN 1993-6 Annex A interaction consume it unchanged.

   Mesh convergence: solve with nSub and 2 nSub subdivisions, report the fine
   solution with meshError = the largest change of max|phi|, max|phi'| and
   max|B| between the two meshes, each normalised by max(its fine-mesh peak,
   1e-3 x its physical scale from the peak torque) so that a quantity that is
   zero in the exact solution cannot read as a mesh error (19 Sep 2026 review
   finding F-B); the caller blocks PASS above TORSION_FE_MESH_BLOCK.

   Units: L, x in mm; E I_w in N.mm4; G I_T in N.mm2; point torque P in N.mm;
   distributed torque w1, w2 in N.mm/mm (= N). Pure functions, no DOM.
   =========================================================================== */
const TORSION_FE_MESH_BLOCK = 0.005;     // 0.5 %: refuse to certify above this
const TORSION_FE_NSUB = 120;             // base subdivisions (doubled once)

function torsionFeKBend(EIw, Le){
  const c=EIw/(Le*Le*Le), L=Le, L2=Le*Le;
  return [[12*c, 6*L*c, -12*c, 6*L*c],[6*L*c, 4*L2*c, -6*L*c, 2*L2*c],
          [-12*c, -6*L*c, 12*c, -6*L*c],[6*L*c, 2*L2*c, -6*L*c, 4*L2*c]];
}
function torsionFeKTors(GIt, Le){
  const c=GIt/(30*Le), L=Le, L2=Le*Le;
  return [[36*c, 3*L*c, -36*c, 3*L*c],[3*L*c, 4*L2*c, -3*L*c, -L2*c],
          [-36*c, -3*L*c, 36*c, -3*L*c],[3*L*c, -L2*c, -3*L*c, 4*L2*c]];
}
/* Consistent nodal loads of a distributed torque varying linearly ta -> tb
   over an element of length Le (Hermite cubic shape functions). */
function torsionFeLoadVector(ta, tb, Le){
  return [Le*(7*ta+3*tb)/20, Le*Le*(3*ta+2*tb)/60, Le*(3*ta+7*tb)/20, -Le*Le*(2*ta+3*tb)/60];
}
/* Banded symmetric positive-definite solve (Cholesky) with half-bandwidth hb.
   K is stored as K[i][j - i + hb] for |j - i| <= hb. Diagonal scaling removes
   the phi / phi' unit disparity exactly as linsolve() does for the beam. */
function torsionFeBandSolve(K, F, hb){
  const n=F.length;
  const scale=new Float64Array(n);
  for(let i=0;i<n;i++){ const d=K[i][hb]; if(!(Number.isFinite(d)&&d>0)) throw new Error('Warping-torsion FE: singular stiffness (no twist restraint?)'); scale[i]=Math.sqrt(d); }
  const C=Array.from({length:n},()=>new Float64Array(hb+1));   // C[i][j-i+hb], j<=i
  for(let i=0;i<n;i++){
    const j0=Math.max(0,i-hb);
    for(let j=j0;j<=i;j++){
      let v=K[i][j-i+hb]/(scale[i]*scale[j]);
      const k0=Math.max(j0,j-hb);
      for(let k=k0;k<j;k++) v-=C[i][k-i+hb]*C[j][k-j+hb];
      if(i===j){ if(!(Number.isFinite(v)&&v>1e-14)) throw new Error('Warping-torsion FE: stiffness matrix is not positive definite (mechanism or ill-conditioned mesh)'); C[i][hb]=Math.sqrt(v); }
      else C[i][j-i+hb]=v/C[j][hb];
    }
  }
  const y=new Float64Array(n), z=new Float64Array(n);
  for(let i=0;i<n;i++){ let v=F[i]/scale[i]; for(let j=Math.max(0,i-hb);j<i;j++) v-=C[i][j-i+hb]*y[j]; y[i]=v/C[i][hb]; }
  for(let i=n-1;i>=0;i--){ let v=y[i]; for(let j=i+1;j<=Math.min(n-1,i+hb);j++) v-=C[j][i-j+hb]*z[j]; z[i]=v/C[i][hb]; }
  return Array.from(z,(v,i)=>v/scale[i]);
}
/* Twist-restraint list of the torsion model: every entry prevents twist
   (phi = 0); warpFix imposes phi' = 0 there as well. Positions in mm. Nothing
   is inferred from the count: a single entry is a torsion cantilever whose
   root warping is fixed only when its flag says so. */
function torsionFeSupports(supports){
  return supports.map(s=>({pos:+s.pos, warpFix:!!s.warpFix}));
}
/* Printed description of the boundary conditions (pure). */
function torsionFeBcText(supports, L){
  const m=v=>(v/1000).toFixed(2).replace(/\.?0+$/,'');
  const sp=torsionFeSupports(supports).slice().sort((p,q)=>p.pos-q.pos);
  const items=sp.map(s=>'x = '+m(s.pos)+' m ('+(s.warpFix? '&phi; = 0, &phi;&prime; = 0: warping fixed' : 'fork, &phi; = 0, warping free')+')');
  const ends=[];
  if(sp[0].pos>1e-6) ends.push('free end at x = 0 (T = B = 0)');
  if(L-sp[sp.length-1].pos>1e-6) ends.push('free end at x = '+m(L)+' m (T = B = 0)');
  return (sp.length===1? 'torsion cantilever: twist restrained at ' : 'twist restrained at ')+items.join(', ')+(ends.length? '; '+ends.join(', ') : '');
}
/* One solve on a mesh of nSub base subdivisions. Returns nodal DOFs and the
   element end forces; see torsionFeRecover for the station values. */
function torsionFeSolveOnce(L, EIw, GIt, supports, torques, nSub){
  const sp=torsionFeSupports(supports);
  const nodes=buildNodes(L, sp, torques, nSub);
  const nN=nodes.length, nd=2*nN, hb=3;
  const idx=new Map(nodes.map((x,i)=>[+x.toFixed(6),i]));
  const K=Array.from({length:nd},()=>new Float64Array(2*hb+1));
  const F=new Float64Array(nd);
  const add=(i,j,v)=>{ K[i][j-i+hb]+=v; };
  const feq=[];   // per element consistent load vector (for the end-force recovery)
  for(let e=0;e<nN-1;e++){
    const Le=nodes[e+1]-nodes[e];
    const kb=torsionFeKBend(EIw,Le), kt=torsionFeKTors(GIt,Le);
    const d=[2*e,2*e+1,2*e+2,2*e+3];
    for(let a=0;a<4;a++) for(let b=0;b<4;b++) add(d[a],d[b],kb[a][b]+kt[a][b]);
    feq.push([0,0,0,0]);
  }
  const nodalP=new Float64Array(nN);
  torques.forEach(ld=>{
    if(ld.type==='point'){ const i=idx.get(+(+ld.pos).toFixed(6)); if(i!=null){ F[2*i]+=ld.P; nodalP[i]+=ld.P; } }
    else if(ld.type==='udl'){
      for(let e=0;e<nN-1;e++){ const xa=nodes[e],xb=nodes[e+1];
        if(xb<=ld.x1+1e-9||xa>=ld.x2-1e-9) continue; const Le=xb-xa;
        const tv=x=>{ if(ld.x2===ld.x1) return ld.w1; const s=(x-ld.x1)/(ld.x2-ld.x1); return ld.w1+(ld.w2-ld.w1)*s; };
        const fv=torsionFeLoadVector(tv(xa),tv(xb),Le);
        for(let a=0;a<4;a++){ F[2*e+a]+=fv[a]; feq[e][a]+=fv[a]; }
      }
    }
  });
  // unconstrained copy for the reactions
  const K0=K.map(r=>Float64Array.from(r)), F0=Float64Array.from(F);
  const fixed=[];
  sp.forEach(s=>{ const i=idx.get(+(+s.pos).toFixed(6)); if(i==null) return; fixed.push(2*i); if(s.warpFix) fixed.push(2*i+1); });
  if(!fixed.length) throw new Error('Warping-torsion FE: no support prevents twist');
  fixed.forEach(r=>{
    for(let k=0;k<2*hb+1;k++){ const c=r+k-hb; if(c<0||c>=nd) continue; K[r][k]=0; K[c][r-c+hb]=0; }
    K[r][hb]=1; F[r]=0;
  });
  const d=torsionFeBandSolve(K,F,hb);
  // reactions R = K0 d - F0 at the constrained DOFs (torque; bimoment when warping is fixed)
  const R=new Float64Array(nd);
  fixed.forEach(r=>{ let s=0; for(let k=0;k<2*hb+1;k++){ const c=r+k-hb; if(c>=0&&c<nd) s+=K0[r][k]*d[c]; } R[r]=s-F0[r]; });
  // element end forces f = K_e d_e - f_eq: [-T(x1), -B(x1), T(x2), B(x2)]
  const ef=[];
  for(let e=0;e<nN-1;e++){
    const Le=nodes[e+1]-nodes[e];
    const kb=torsionFeKBend(EIw,Le), kt=torsionFeKTors(GIt,Le);
    const de=[d[2*e],d[2*e+1],d[2*e+2],d[2*e+3]];
    const f=[0,0,0,0];
    for(let a=0;a<4;a++){ let s=0; for(let b=0;b<4;b++) s+=(kb[a][b]+kt[a][b])*de[b]; f[a]=s-feq[e][a]; }
    ef.push({TL:-f[0], BL:-f[1], TR:f[2], BR:f[3]});
  }
  const reactions=sp.map(s=>{ const i=idx.get(+(+s.pos).toFixed(6)); return {pos:s.pos, warpFix:s.warpFix, T:i!=null? R[2*i]:0, B:(i!=null&&s.warpFix)? R[2*i+1]:0}; });
  const nodeOf=s=>idx.get(+(+s.pos).toFixed(6));
  return {nodes,d,ef,nodalP,reactions,nElem:nN-1,
    supportNodes:new Set(sp.map(nodeOf).filter(i=>i!=null)),
    warpNodes:new Set(sp.filter(s=>s.warpFix).map(nodeOf).filter(i=>i!=null))};   // phi' fixed: B jumps by the reaction bimoment
}
/* Station values from one solve (see the header for the recovery rules). */
function torsionFeRecover(sol, EIw, GIt){
  const {nodes,d,ef,nodalP,supportNodes}=sol, warpNodes=sol.warpNodes||new Set();
  const nN=nodes.length;
  const xs=[],phi=[],p1=[],p2=[],p3=[],T=[],B=[];
  const push=(x,ph,dph,Bv,Tv)=>{ xs.push(x); phi.push(ph); p1.push(dph); B.push(Bv); T.push(Tv); p2.push(Bv/EIw); p3.push((GIt*dph-Tv)/EIw); };
  for(let i=0;i<nN;i++){
    const x=nodes[i], ph=d[2*i], dph=d[2*i+1];
    const left = i>0? ef[i-1] : null, right = i<nN-1? ef[i] : null;
    // B is continuous at an ordinary interior node (average the two element
    // values); at an interior warping-fixed node the reaction bimoment makes
    // B jump, so the one-sided element values are kept (larger magnitude at x)
    const bJump = !!(left&&right && warpNodes.has(i));
    const BLv = left? left.BR : null, BRv = right? right.BL : null;
    const Bv = bJump? (Math.abs(BLv)>=Math.abs(BRv)? BLv : BRv) : (left&&right? 0.5*(BLv+BRv) : (left? BLv : BRv));
    const TLv = left? left.TR : null, TRv = right? right.TL : null;
    const jump = (left&&right && (Math.abs(nodalP[i])>1e-9 || supportNodes.has(i))) || bJump;
    if(jump){
      if(x-0.01>xs[xs.length-1]) push(x-0.01, ph, dph, bJump? BLv : Bv, TLv);
      push(x, ph, dph, Bv, Math.abs(TLv)>=Math.abs(TRv)? TLv : TRv);
      push(x+0.01, ph, dph, bJump? BRv : Bv, TRv);
    } else push(x, ph, dph, Bv, TLv!=null? TLv : TRv);
  }
  return {xs,phi,p1,p2,p3,T,B};
}
function torsionFeMax(arr){ let m=0; arr.forEach(v=>{ if(Math.abs(v)>m) m=Math.abs(v); }); return m; }
/* Public entry. opts = {L, EIw, GIt, supports:[{pos (mm), warpFix?}], torques:
   [{type:'point', pos, P} | {type:'udl', x1, x2, w1, w2}], nSub?, refine?}.
   Returns the fine-mesh solution in p385Solve's shape plus T (total torque),
   B (bimoment), reactions, nElem, nElemCoarse, meshError, meshErrorParts,
   converged, X = L/a, bc (support list), method label. */
function warpingTorsionFE(opts){
  const {L,EIw,GIt,supports,torques}=opts;
  if(!(L>0)) throw new Error('Warping-torsion FE: L must be > 0');
  if(!(EIw>0)) throw new Error('Warping-torsion FE: E I_w must be > 0');
  if(!(GIt>0)) throw new Error('Warping-torsion FE: G I_T must be > 0');
  const nSub=opts.nSub||TORSION_FE_NSUB;
  const refine=opts.refine!==false;
  const aa=Math.sqrt(EIw/GIt);
  const coarse=torsionFeSolveOnce(L,EIw,GIt,supports,torques,nSub);
  const rc=torsionFeRecover(coarse,EIw,GIt);
  let fine=coarse, rf=rc, meshError=0, parts={phi:0,p1:0,B:0};
  if(refine){
    fine=torsionFeSolveOnce(L,EIw,GIt,supports,torques,2*nSub);
    rf=torsionFeRecover(fine,EIw,GIt);
    // Mesh measure: the change of each peak on doubling the mesh, normalised
    // by max(the fine-mesh peak, 1e-3 x the physical scale of that quantity
    // set by the peak torque T_max): T_max L/GI_T (St Venant twist), T_max/GI_T
    // (rate of twist), T_max min(a, L) (bimoment: T a/2 for a long member, T L/4
    // for a short one). A peak below a thousandth of its scale is identically
    // zero in the exact solution up to round-off (the bimoment of a warping-
    // free cantilever under a tip torque, pure St Venant), so its round-off on
    // the two meshes no longer reads as a spurious "mesh error" (19 Sep 2026
    // review finding F-B); a physically significant peak keeps the plain
    // relative measure.
    const Tmax=torsionFeMax(rf.T);
    const rel=(a,b,scale)=>{ const m=torsionFeMax(b), s=Math.max(m,scale); return s>1e-300? Math.abs(torsionFeMax(a)-m)/s : 0; };
    parts={phi:rel(rc.phi,rf.phi,1e-3*Tmax*L/GIt), p1:rel(rc.p1,rf.p1,1e-3*Tmax/GIt), B:rel(rc.B,rf.B,1e-3*Tmax*Math.min(aa,L))};
    meshError=Math.max(parts.phi,parts.p1,parts.B);
  }
  return Object.assign({},rf,{X:L/aa, aa, nElem:fine.nElem, nElemCoarse:coarse.nElem, meshError, meshErrorParts:parts,
    converged:meshError<=TORSION_FE_MESH_BLOCK, reactions:fine.reactions, bc:torsionFeSupports(supports),
    method:'fe', methodLabel:'warping-torsion FE ('+fine.nElem+' elements)'});
}


/* ==== beam-v03 module 13 ==== */
/* ===========================================================================
   5. DIAGRAMS (inline SVG)
   =========================================================================== */
function svgEl(W,H,inner,attrs){ return `<svg class="diag${attrs?' diag-hover':''}" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg"${attrs||''}>${inner}</svg>`; }
function plot(xs,ys,opt){
  const W=540,H=150,padL=46,padR=16,padT=16,padB=24;
  const xmax=xs[xs.length-1]||1;
  const ymax=Math.max(1e-9,...ys.map(v=>Math.abs(v)));
  const X=x=>padL+x/xmax*(W-padL-padR);
  const midY=padT+(H-padT-padB)/2;
  const amp=(H-padT-padB)/2-4;
  const sgn=opt.flip?-1:1;
  const Y=v=>midY - sgn*v/ymax*amp;
  let path=`M ${X(xs[0])} ${Y(ys[0])}`;
  for(let i=1;i<xs.length;i++) path+=` L ${X(xs[i])} ${Y(ys[i])}`;
  let fill=`M ${X(xs[0])} ${midY} L ${X(xs[0])} ${Y(ys[0])}`;
  for(let i=1;i<xs.length;i++) fill+=` L ${X(xs[i])} ${Y(ys[i])}`;
  fill+=` L ${X(xs[xs.length-1])} ${midY} Z`;
  // peak
  let pk=0,pi=0; ys.forEach((v,i)=>{ if(Math.abs(v)>Math.abs(pk)){pk=v;pi=i;} });
  const peakLbl=`${opt.fmt(pk)} ${opt.unit}`;
  const lx=X(xs[pi]), ly=Y(pk);
  let inner=`
    <line x1="${padL}" y1="${midY}" x2="${W-padR}" y2="${midY}" stroke="#888" stroke-width="1"/>
    <path d="${fill}" fill="${opt.fill}" opacity="0.5"/>
    <path d="${path}" fill="none" stroke="${opt.color}" stroke-width="2"/>
    <circle cx="${lx}" cy="${ly}" r="3" fill="${opt.color}"/>
    <text x="${Math.min(Math.max(lx,padL+30),W-padR-30)}" y="${ly+ (sgn*pk>=0? -6:14)}" font-family="Arial" font-size="11" font-weight="700" fill="${opt.color}" text-anchor="middle">${peakLbl}</text>
    <text x="6" y="${midY+4}" font-family="Arial" font-size="10" fill="#555">0</text>
    <text x="${padL}" y="${H-8}" font-family="Arial" font-size="10" fill="#555">0</text>
    <text x="${W-padR}" y="${H-8}" font-family="Arial" font-size="10" fill="#555" text-anchor="end">${xmax.toFixed(2)} m</text>`;
  // 20 Sep 2026 hover readout (owner: "see values on the diagram at any point
  // when I hover"). The sample arrays (7 s.f., x in m; 20 Sep 2026 review: 4
  // s.f. let the 2-dp readout disagree with the peak label and the forces
  // table above ~1000 kN.m, e.g. 2607.00 for 2606.92) and every number of the
  // data -> viewBox mapping above ride on the <svg> as data-* attributes, and a
  // hidden <g class="hover"> (guide line, marker, haloed readout) waits at the
  // end so it paints on top. installDiagramHover() below drives it through ONE
  // delegated listener set - no per-svg state, so innerHTML re-renders are free.
  // opt.name ('V','M','delta','T') and opt.caption feed the readout only; with
  // no name the readout is "x = .. m   <value> <unit>" except for the three
  // unambiguous defaults: unit 'kN' -> V, unit 'mm' -> delta, and a flipped
  // (sagging-down) moment plot -> M (torsion is not flipped, so it stays blank).
  const sig=v=>+Number(v).toPrecision(7), esc=v=>String(v==null?'':v).replace(/"/g,'&quot;').replace(/</g,'&lt;');
  const unit=opt.unit||'';
  const name=opt.name!=null? String(opt.name) : unit==='kN'?'V' : unit==='mm'?'δ' : (opt.flip&&/^kN.{0,9}m$/.test(unit))?'M':'';
  const attrs=` data-xs="${JSON.stringify(xs.map(sig))}" data-ys="${JSON.stringify(ys.map(sig))}"`+
    ` data-unit="${esc(unit)}" data-fmt-dp="${opt.dp!=null?+opt.dp:2}" data-xlabel="${esc(opt.xlabel||'x')}" data-name="${esc(name)}" data-caption="${esc(opt.caption||'')}"`+
    ` data-padl="${padL}" data-padr="${padR}" data-padt="${padT}" data-padb="${padB}" data-xmax="${sig(xmax)}" data-ymax="${sig(ymax)}" data-flip="${opt.flip?1:0}" data-midy="${midY}" data-amp="${amp}"`;
  inner+=`
    <g class="hover" style="display:none" pointer-events="none">
      <line class="hover-x" x1="${padL}" y1="${padT}" x2="${padL}" y2="${H-padB}" stroke="#444" stroke-width="1" stroke-dasharray="3 2"/>
      <circle class="hover-pt" cx="${padL}" cy="${midY}" r="4" fill="${opt.color}" stroke="#fffdf8" stroke-width="1.5"/>
      <text class="hover-txt" x="${padL+8}" y="${midY-9}" font-family="Arial" font-size="11" font-weight="700" fill="#111" stroke="#fffdf8" stroke-width="4" paint-order="stroke" stroke-linejoin="round" text-anchor="start"></text>
    </g>`;
  return svgEl(W,H,inner,attrs);
}

/* ---- hover readout driver (20 Sep 2026) ------------------------------------
   Pure helpers first (they run in the vm test harness, which has no DOM):
   diagHoverData(svg)      -> the data-* attributes parsed once per element
                              (cached in a WeakMap keyed by the svg, so a
                              re-rendered svg simply parses afresh);
   diagHoverReadout(d,vx)  -> for a viewBox x: the interpolated sample, the
                              marker position on the curve, the readout text
                              and its anchor/position kept inside the viewBox.
   installDiagramHover(root) (default document) attaches ONE delegated listener
   set - mousemove / mouseleave (capture, it does not bubble) / touchmove /
   touchend / touchcancel - for every svg.diag-hover under root, present now or
   rendered later by innerHTML. Idempotent per root; returns false when the
   root cannot take listeners (the harness stub document). No inline handlers,
   no eval: CSP-safe. */
function diagHoverData(svg){
  const cache=diagHoverData.cache||(diagHoverData.cache=(typeof WeakMap==='function')? new WeakMap() : null);
  if(cache&&cache.has(svg)) return cache.get(svg);
  const g=k=>svg.getAttribute('data-'+k), n=k=>+g(k);
  let xs,ys; try{ xs=JSON.parse(g('xs')||'[]'); ys=JSON.parse(g('ys')||'[]'); }catch(e){ xs=[]; ys=[]; }
  const vb=String(svg.getAttribute('viewBox')||'0 0 540 150').trim().split(/[\s,]+/).map(Number);
  const d={xs,ys,W:vb[2]||540,H:vb[3]||150,padL:n('padl'),padR:n('padr'),padT:n('padt'),padB:n('padb'),xmax:n('xmax')||1,ymax:n('ymax')||1e-9,
    flip:g('flip')==='1',midY:n('midy'),amp:n('amp'),unit:g('unit')||'',name:g('name')||'',dp:(g('fmt-dp')!=null&&g('fmt-dp')!=='')? +g('fmt-dp') : 2,
    xlabel:g('xlabel')||'x',caption:g('caption')||''};
  if(!(d.midY>0)) d.midY=d.padT+(d.H-d.padT-d.padB)/2;
  if(!(d.amp>0)) d.amp=(d.H-d.padT-d.padB)/2-4;
  if(cache) cache.set(svg,d);
  return d;
}
function diagHoverReadout(d,vx){
  const n=d.xs.length; if(!n||d.ys.length!==n) return null;
  const span=d.W-d.padL-d.padR;
  // viewBox x -> member x, clamped to the sampled range
  let x=(vx-d.padL)/span*d.xmax; x=Math.min(Math.max(x,d.xs[0]),d.xs[n-1]);
  // segment [lo,hi] by bisection, linear interpolation inside it (a zero-length
  // segment - the two samples either side of a jump - takes its first end)
  let lo=0,hi=n-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(d.xs[m]<=x) lo=m; else hi=m; }
  const x0=d.xs[lo],x1=d.xs[hi]; let t=(x1>x0)? (x-x0)/(x1-x0) : 0; t=Math.min(Math.max(t,0),1);
  const y=d.ys[lo]+(d.ys[hi]-d.ys[lo])*t;
  const cx=d.padL+x/d.xmax*span, cy=d.midY-(d.flip?-1:1)*y/d.ymax*d.amp;
  let ys=y.toFixed(d.dp); if(/^-0\.?0*$/.test(ys)) ys=ys.slice(1);   // no "-0.00"
  const text=(d.caption? d.caption+': ' : '')+`${d.xlabel} = ${x.toFixed(2)} m ${d.name? d.name+' = ' : ''}${ys} ${d.unit}`.trim();
  // keep the readout inside the viewBox: flip the anchor near the right edge,
  // and sit the text on the EMPTY side of the axis at this x (the curve is on
  // one side only there), tied to the marker by the guide line - so it never
  // fights the peak label, which lives on the curve's outer side
  const est=text.length*6.3;
  const flipAnchor=cx+8+est>d.W-2;
  const tx=flipAnchor? cx-8 : cx+8, anchor=flipAnchor? 'end' : 'start';
  const ty=(cy>=d.midY)? d.midY-12 : d.midY+20;
  return {x,y,cx,cy,text,tx,ty,anchor};
}
function diagHoverPoint(svg,clientX,clientY){
  // pointer -> viewBox units: the screen CTM when the browser gives one, else
  // the bounding-rect ratio (width:100%, default xMidYMid meet => uniform scale)
  try{
    if(typeof svg.getScreenCTM==='function'){
      const m=svg.getScreenCTM();
      if(m){ const inv=m.inverse();
        if(typeof DOMPoint==='function'){ const p=new DOMPoint(clientX,clientY).matrixTransform(inv); return {x:p.x,y:p.y}; }
        if(typeof svg.createSVGPoint==='function'){ let p=svg.createSVGPoint(); p.x=clientX; p.y=clientY; p=p.matrixTransform(inv); return {x:p.x,y:p.y}; }
      }
    }
  }catch(e){ /* fall through to the ratio */ }
  const d=diagHoverData(svg), r=svg.getBoundingClientRect();
  const sc=Math.min(r.width/d.W,r.height/d.H)||1;
  return {x:(clientX-r.left-(r.width-d.W*sc)/2)/sc, y:(clientY-r.top-(r.height-d.H*sc)/2)/sc};
}
function diagHoverShow(svg,clientX,clientY){
  const grp=svg.querySelector('.hover'); if(!grp) return;
  const d=diagHoverData(svg), p=diagHoverPoint(svg,clientX,clientY), r=diagHoverReadout(d,p.x);
  if(!r){ grp.style.display='none'; return; }
  const ln=grp.querySelector('.hover-x'), pt=grp.querySelector('.hover-pt'), tx=grp.querySelector('.hover-txt');
  if(ln){ ln.setAttribute('x1',r.cx); ln.setAttribute('x2',r.cx); }
  if(pt){ pt.setAttribute('cx',r.cx); pt.setAttribute('cy',r.cy); }
  if(tx){ tx.setAttribute('x',r.tx); tx.setAttribute('y',r.ty); tx.setAttribute('text-anchor',r.anchor); tx.textContent=r.text; }
  grp.style.display='';
}
function diagHoverHide(svg){ const grp=svg&&svg.querySelector&&svg.querySelector('.hover'); if(grp) grp.style.display='none'; }
function installDiagramHover(root){
  root=root||(typeof document!=='undefined'? document : null);
  if(!root||typeof root.addEventListener!=='function') return false;
  const roots=installDiagramHover.roots||(installDiagramHover.roots=[]);
  if(roots.indexOf(root)>=0) return false;   // idempotent: one listener set per root
  roots.push(root);
  let active=null;   // the svg currently showing a readout (dropped when detached)
  const svgOf=t=>(t&&typeof t.closest==='function')? t.closest('svg.diag-hover') : null;
  const hideActive=()=>{ if(active){ diagHoverHide(active); active=null; } };
  const move=(target,clientX,clientY)=>{
    const svg=svgOf(target);
    if(!svg){ hideActive(); return; }
    if(active&&active!==svg) hideActive();
    active=svg; diagHoverShow(svg,clientX,clientY);
  };
  root.addEventListener('mousemove',e=>move(e.target,e.clientX,e.clientY),{passive:true});
  root.addEventListener('mouseleave',e=>{ const svg=svgOf(e.target); if(svg){ diagHoverHide(svg); if(active===svg) active=null; } },{capture:true,passive:true});
  root.addEventListener('touchmove',e=>{ const t=e.touches&&e.touches[0]; if(t) move(t.target||e.target,t.clientX,t.clientY); },{passive:true});
  root.addEventListener('touchend',hideActive,{passive:true});
  root.addEventListener('touchcancel',hideActive,{passive:true});
  return true;
}
function beamDiagram(a){
  // Collision-free loading sketch:
  //  - distributed loads whose extents overlap are STACKED in tiers, each in
  //    its own colour with its label on its own tier line;
  //  - point/moment load labels live in a collision-managed band above the
  //    load zone (two rows, side-nudged like the section view);
  //  - all text carries a paper halo, and the canvas height grows with the
  //    number of tiers so nothing is ever forced on top of anything else.
  const W=540,padL=46,padR=16;
  const L=S.L, xmax=L||1;
  const X=x=>padL+x/xmax*(W-padL-padR);
  const halo=' stroke="#fffdf8" stroke-width="3" paint-order="stroke" stroke-linejoin="round"';
  const distCols=['#06c','#0a7f5a','#c2620a','#8b3fc9','#0a7fa8','#b00'];
  const loads=S.loads;
  // ---- tier assignment for distributed loads (overlap => next tier up) ----
  const dists=[]; loads.forEach((ld,i)=>{ if(ld.type==='udl'||ld.type==='trap') dists.push({ld,i}); });
  const tierRanges=[];
  dists.forEach(o=>{
    const x1=Math.min(+o.ld.x1,+o.ld.x2), x2=Math.max(+o.ld.x1,+o.ld.x2);
    let t=0;
    while(tierRanges[t] && tierRanges[t].some(r=>x1<r[1]-1e-9 && r[0]<x2-1e-9)) t++;
    (tierRanges[t]=tierRanges[t]||[]).push([x1,x2]);
    o.tier=t; o.x1=x1; o.x2=x2;
  });
  const nT=Math.max(tierRanges.length,0), tierH=24;
  const havePts=loads.some(ld=>ld.type==='point'||ld.type==='moment');
  const topMargin=6;
  const bandRows=havePts?2:0, bandRowY=[topMargin+10,topMargin+23];
  const distTop=topMargin+(havePts?bandRows*13+4:2);
  const yB=distTop + (nT>0? nT*tierH : 34) + 8;
  const H=yB+72;
  let inner=`<line x1="${X(0)}" y1="${yB}" x2="${X(L)}" y2="${yB}" stroke="#111" stroke-width="3"/>`;
  // end supports from the in-plane type of each end: pinned = triangle, fixed =
  // hatched wall, guided (rotation held, vertical free) = sliding block, free = nothing
  endsToSupports(S).forEach(sp=>{ const x=X(sp.pos), sgn=sp.end===1? -1 : 1;
    if(sp.type==='pinned'){
      inner+=`<polygon points="${x},${yB} ${x-7},${yB+13} ${x+7},${yB+13}" fill="none" stroke="#111" stroke-width="1.6"/>
        <line x1="${x-10}" y1="${yB+16}" x2="${x+10}" y2="${yB+16}" stroke="#111" stroke-width="1.4"/>`;
    } else if(sp.type==='fixed'){
      inner+=`<line x1="${x}" y1="${yB-16}" x2="${x}" y2="${yB+16}" stroke="#111" stroke-width="2.4"/>`;
      for(let k=-14;k<=14;k+=6) inner+=`<line x1="${x}" y1="${yB+k}" x2="${x+sgn*7}" y2="${yB+k+6}" stroke="#111" stroke-width="1"/>`;
    } else if(sp.type==='guided'){
      inner+=`<rect x="${x-3}" y="${yB-14}" width="6" height="28" fill="none" stroke="#111" stroke-width="1.6"/>
        <line x1="${x+sgn*7}" y1="${yB-16}" x2="${x+sgn*7}" y2="${yB+16}" stroke="#111" stroke-width="2.4"/>`;
      for(let k=-14;k<=14;k+=6) inner+=`<line x1="${x+sgn*7}" y1="${yB+k}" x2="${x+sgn*14}" y2="${yB+k+6}" stroke="#111" stroke-width="1"/>`;
    }
  });
  // ---- distributed loads, tier by tier ----
  dists.forEach((o,di)=>{
    const ld=o.ld, col=distCols[di%distCols.length];
    const x1=X(o.x1), x2=X(o.x2);
    const top=yB-8-(o.tier+1)*tierH+tierH-16;   // block top for this tier
    inner+=`<line x1="${x1}" y1="${top}" x2="${x2}" y2="${top}" stroke="${col}" stroke-width="1.6"/>`;
    const N=Math.max(2,Math.round((x2-x1)/22));
    for(let k=0;k<=N;k++){ const xx=x1+(x2-x1)*k/N;
      inner+=`<line x1="${xx}" y1="${top}" x2="${xx}" y2="${yB-2}" stroke="${col}" stroke-width="1.1" opacity="${o.tier>0?0.55:0.9}"/>
        <polygon points="${xx},${yB-2} ${xx-3},${yB-9} ${xx+3},${yB-9}" fill="${col}" opacity="${o.tier>0?0.6:1}"/>`; }
    const w1=ld.type==='trap'?ld.w1:ld.w, w2=ld.type==='trap'?ld.w2:ld.w;
    const lbl=(ld.type==='trap'?`${w1} to ${w2}`:`${ld.w}`)+' kN/m'+(ld.case?` (${ld.case})`:'');
    const lx=Math.min(Math.max((x1+x2)/2,padL+34),W-padR-34);
    inner+=`<text x="${lx}" y="${top-4}" font-family="Arial" font-size="10.5" font-weight="700" fill="${col}" text-anchor="middle"${halo}>${lbl}</text>`;
  });
  // ---- point + moment loads: arrows to the beam, labels in the top band ----
  const bandItems=[];
  loads.forEach(ld=>{
    if(ld.type==='point'){ const x=X(+ld.pos); const dir=(ld.P>=0)?1:-1;
      if(dir>0){
        const startY=distTop+2;
        inner+=`<line x1="${x}" y1="${startY}" x2="${x}" y2="${yB-2}" stroke="#b00" stroke-width="2"/>
          <polygon points="${x},${yB-2} ${x-4},${yB-12} ${x+4},${yB-12}" fill="#b00"/>`;
        bandItems.push({x, text:`${Math.abs(ld.P)} kN`, color:'#b00', anchorY:startY});
      } else { // uplift: below the beam
        inner+=`<line x1="${x}" y1="${yB+34}" x2="${x}" y2="${yB+2}" stroke="#b00" stroke-width="2"/>
          <polygon points="${x},${yB+2} ${x-4},${yB+12} ${x+4},${yB+12}" fill="#b00"/>
          <text x="${x}" y="${yB+46}" font-family="Arial" font-size="10.5" font-weight="700" fill="#b00" text-anchor="middle"${halo}>${Math.abs(ld.P)} kN</text>`;
      }
    } else if(ld.type==='moment'){ const x=X(+ld.pos);
      inner+=`<path d="M ${x-10} ${yB-22} A 11 11 0 1 1 ${x-11} ${yB-20}" fill="none" stroke="#7a4" stroke-width="2"/>
        <polygon points="${x-11},${yB-20} ${x-15},${yB-24} ${x-7},${yB-26}" fill="#7a4"/>`;
      bandItems.push({x, text:`${Math.abs(ld.M)} kN·m`, color:'#5a7a2a', anchorY:yB-26});
    }
  });
  if(bandItems.length){
    const est=t=>t.length*6.0+6;
    bandItems.sort((p,q)=>p.x-q.x);
    const lastEnd=[-1e9,-1e9];
    bandItems.forEach(it=>{
      const w=est(it.text);
      let row=0, lab=it.x;
      if(it.x-w/2 <= lastEnd[0]+4){
        if(it.x-w/2 > lastEnd[1]+4) row=1;
        else { row=lastEnd[0]<=lastEnd[1]?0:1; lab=lastEnd[row]+4+w/2; }
      }
      lab=Math.min(Math.max(lab,w/2+2),W-w/2-2);
      lastEnd[row]=lab+w/2;
      const y=bandRowY[row];
      if(Math.abs(lab-it.x)>6) inner+=`<line x1="${lab}" y1="${y+2}" x2="${it.x}" y2="${it.anchorY}" stroke="${it.color}" stroke-width="0.7" opacity="0.7"/>`;
      inner+=`<text x="${lab}" y="${y}" font-family="Arial" font-size="10.5" font-weight="700" fill="${it.color}" text-anchor="middle"${halo}>${it.text}</text>`;
    });
  }
  inner+=`<line x1="${X(0)}" y1="${H-14}" x2="${X(L)}" y2="${H-14}" stroke="#999" stroke-width="1"/>
    <text x="${(X(0)+X(L))/2}" y="${H-4}" font-family="Arial" font-size="10.5" fill="#555" text-anchor="middle">L = ${L.toFixed(2)} m</text>`;
  return svgEl(W,H,inner);
}


/* ==== beam-v03 module 14 ==== */
/* ===========================================================================
   6. REPORT RENDER
   render(): verdict banner, then (EC3, since 20 Sep 2026) the reference software
   brief of js/06-brief-reference.js alone, or (BS 5950) the long CED
   report built here.
   =========================================================================== */
const f1=(v,d=1)=>{ if(!isFinite(v))return" "; const n=Math.abs(v)<5e-7?0:v; return n.toFixed(d); };
  const g=(v,d=2)=>{ if(!isFinite(v))return" "; const s=(Math.abs(v)<5e-7?0:v).toFixed(d); return s.includes('.') ? s.replace(/(\.\d*?)0+$/,'$1').replace(/\.$/,'') : s; };
const sname=k=>sectionDisplayName(k);   // "200 x 75 x 23" for every family (js/03-state-ui.js; 20 Sep 2026: the old /x/g -> ' ' printed "200 75 23")
function st(ok,t,fl){ return `<div class="status ${ok?'ok':'fail'}">${ok?t:(fl||'FAIL')}</div>`; }

function render(){
  const rep=$("report");
  let a,c;
  try{ a=analyse(); c=checks(a); }
  catch(err){ rep.innerHTML=`<div class="err">Could not analyse: ${err}</div>`; return; }
  const sec=a.sec;
  const sci = S.code==='EC3' && (S.restraint||'full')==='full';
  const sciU = S.code==='EC3' && !sci;
  // EC3 unrestrained Mcr method actually used by the check engine (c.mcrMethod is
  // set by both routes; 'standard' = closed form, 'eigen' = FE eigensolver)
  const mcrStd = sciU && c.mcrMethod==='standard';
  const famLabel = sec.isBox? (S.family==='rhs'? `RHS [${sec.boxType==='CF'?'Cold-formed':'Hot-finished'}]` : `SHS [${sec.boxType==='CF'?'Cold-formed':'Hot-finished'}]`) : S.family==='ub'? 'UB' : S.family==='uc'? 'UC' : 'PFC';
  const gradeTxt=`${S.grade} (p<sub>y</sub> = ${g(a.py,0)} N/mm )`;
  // load summary lines
  const swE=selfWeightEccentricity(sec);
  const swOff=(S.eccOn && Math.abs(swE)>1e-9)?`  [e=${g(swE,1)} mm, z_g=0 mm]`:'';
  const swPl = typeof plateGeom==='function'? plateGeom(sec) : null;
  const autoSwLine=`Self-weight  ${g(selfWeightValue(sec),4)} kN/m ?  0 ${g(S.L)} m  [Dead (G)]${swOff}  [automatic${swPl? ', incl. '+g(swPl.w,0)+'&times;'+g(swPl.t,0)+' '+(swPl.side==='top'?'top':'bottom')+' plate':''}]`;
  const userLoadLines=S.loads.filter(ld=>!ld.isSelfWeight).map(ld=>{
    const sw = ld.isSelfWeight? '  [self-weight]' : '';
    const cs = `  [${CASE_LABELS[ld.case]||ld.case}]`;
    const off = (S.eccOn && ld.type!=='moment') ? `  [e=${g(ld.e||0,0)} mm${loadHeightPerLoadOn()?`, z_g=${g(loadZgValue(ld),0)} mm`:''}]` : '';
    if(ld.type==='point') return `Point  ${g(ld.P)} kN ?  @ ${g(ld.pos)} m${cs}${off}${sw}`;
    if(ld.type==='moment') return `Moment ${g(ld.M)} kN m ?  @ ${g(ld.pos)} m${cs}${sw}`;
    if(ld.type==='udl') return `UDL    ${g(ld.w)} kN/m ?  ${g(ld.x1)} ${g(ld.x2)} m${cs}${off}${sw}`;
    if(ld.type==='trap') return `Trap   ${g(ld.w1)}?${g(ld.w2)} kN/m ?  ${g(ld.x1)} ${g(ld.x2)} m${cs}${off}${sw}`;
    return "";
  });
  const loadLines=[autoSwLine,...userLoadLines].join("<br>");
  const sectionView = typeof sectionLoadLineView === 'function' ? sectionLoadLineView(sec) : '';
  // End 1 / End 2 reactions of the governing-moment combination: R (upward positive) and the
  // end moment in the diagram convention as the end type carries them; "guided: M only"
  const reactLine=endsList().map(e=>{
    const r=a.reactions.find(q=>q.end===e.n);
    if(!r) return `End ${e.n} (x = ${g(e.x)} m): free end, no reaction`;
    if(r.type==='guided') return `End ${e.n} (x = ${g(e.x)} m): guided: M only, M = ${f1(reactionEndMomentKNm(r),2)} kN m`;
    return `End ${e.n} (x = ${g(e.x)} m): R = ${f1(r.V/1000,2)} kN`+(r.type==='fixed'? `, M = ${f1(reactionEndMomentKNm(r),2)} kN m` : ' (pinned: no M)');
  }).join(" &nbsp; | &nbsp; ");
  // End conditions line (19 Sep 2026 scope): the seven DOF flags of each end, the derived end types, the preset and the hinges
  const endsLine=`<div class="note" style="margin-left:0">End conditions: ${endsConditionsLine(S)}.</div>`;
  // uplift / hold-down rows (every combination's reactions; item 1.2)
  const HD=c.holdDown||null;
  const upliftLines=(HD&&HD.rows&&HD.rows.length)
    ? HD.rows.map(u=>{
        if(u.level==='sls') return `<div class="note" style="margin-left:0"><b>Uplift at SLS only, End ${u.n} (x = ${g(u.pos/1000,2)} m):</b> R = &minus;${f1(Math.abs(u.RSls),2)} kN (SLS combination ${u.comboSls}); no ULS combination lifts this end, including the &gamma;<sub>G,inf</sub> companions with G at 1.0 (STR set B) and 0.9 (EQU set A) (advisory).</div>`;
        return `<div class="note" style="margin-left:0;color:${u.holdDown?'#374151':'#b91c1c'}"><b>Hold-down ${u.holdDown?'provided':'required'} at End ${u.n} (x = ${g(u.pos/1000,2)} m):</b> R = &minus;${f1(Math.abs(u.RUls),2)} kN (combination ${u.comboUls})${u.RSls!=null? `; SLS uplift &minus;${f1(Math.abs(u.RSls),2)} kN (${u.comboSls})`:''}; ${u.nCombos} combination(s) lift this end${u.holdDown? ' &mdash; design the hold-down connection for this force (advisory)' : ' &mdash; NOT VERIFIED until "hold-down provided" is ticked for this end'}.</div>`; }).join('')
    : (a.uplift? `<div class="note" style="margin-left:0">Uplift: no end lifts in any of the ${a.uplift.nCombos!=null? a.uplift.nCombos : a.ulsResults.length+a.slsResults.length} combinations (all reactions &ge; 0${(a.ulsCompanions&&a.ulsCompanions.length)? '; incl. the '+a.ulsCompanions.length+' &gamma;<sub>G,inf</sub> companions with G at 1.0 and 0.9' : ''}).</div>` : '');
  // gamma_G,inf companions (reactions only)
  const patternBlock=a.companionNote? `<div class="note" style="margin-left:0">${a.companionNote}</div>` : '';

  // load combination results table   one row per enabled combo, flagging which one governs
  const comboRows = a.ulsResults.map(res=>{
    const isVgov = res===a.governV, isMgov = res===a.governM;
    const flags = [isVgov?'governs shear':null, isMgov?'governs moment':null].filter(Boolean).join(', ');
    return `<tr><td>${res.combo.label}</td><td class="num">${f1(res.Vmax/1000,2)}</td><td class="num">${f1(res.Mmax/1e6,2)} @ ${g(res.Mpos/1000,2)}</td><td style="font-weight:${flags?'700':'400'};color:${flags?'#166534':'#666'}">${flags||' '}</td></tr>`;
  }).join('');
  const slsRows = a.slsResults.map(res=>{
    const isDgov = res===a.governD;
    return `<tr><td>${res.combo.label}</td><td class="num">${f1(res.dmax,2)}</td><td style="font-weight:${isDgov?'700':'400'};color:${isDgov?'#166534':'#666'}">${isDgov?'governs deflection':' '}</td></tr>`;
  }).join('');
  const combosBlock = `
  <div class="section-title smallgap">Load Combinations Considered</div>
  ${patternBlock}
  <table class="force-table" style="font-size:13px">
    <thead><tr><th>ULS combination</th><th>Max F<sub>v</sub> (kN)</th><th>Max M<sub>x</sub> (kN m @ m)</th><th>Governs</th></tr></thead>
    <tbody>${comboRows}</tbody>
  </table>
  <table class="force-table" style="font-size:13px;margin-top:4px">
    <thead><tr><th>SLS combination</th><th>Deflection (mm)</th><th>Governs</th></tr></thead>
    <tbody>${slsRows}</tbody>
  </table>`;

  const verdict=c.pass?"PASS":c.utils.some(u=>!Number.isFinite(u.val)||u.val>1.0001)?"FAIL":"NOT VERIFIED";
  const verdictCombo=c.gov.combo || (c.gov.name==='Deflection'?a.governD.combo.label:
    c.gov.name.startsWith('Shear')?a.governV.combo.label:
    c.gov.name.startsWith('Member buckling')&&c.buck?c.buck.combo||a.governM.combo.label:
    c.gov.name.startsWith('LTB ')&&c.ltb?c.ltb.governCombo||'':'');
  const codeLabel = S.code==='EC3'? (sci? 'EN 1993-1-1 (UK NA) &mdash; fully restrained beam' : (c.ltb&&c.ltb.na? 'EN 1993-1-1 (UK NA) &mdash; closed section beam' : 'EN 1993-1-1 (UK NA) &mdash; unrestrained beam (LTB, M<sub>cr</sub> '+(mcrStd?'standard closed form':'FE eigenvalue')+')')) : 'BS 5950-1:2000';
  const banner=`<div class="banner ${c.pass?'pass':'failb'}">
    <div><div class="verdict">${verdict}</div><div style="font-size:12px;color:#374151;font-family:Arial">${codeLabel} member check   ${sname(sec.key)} ${famLabel}   ${S.grade}</div></div>
    <div class="util">Governing: <b>${c.gov.name} = ${g(c.gov.val,3)}</b>${verdictCombo?' ('+verdictCombo+')':''}<br>
      ${sci
        ? `Shear ${g(c.utils[0].val,2)} &bull; Bending ${g(c.utils[1].val,2)} &bull; Defl ${g(c.utils[2].val,2)}`
        : sciU && c.ltb && c.ltb.na
        ? `Shear ${g(c.utils[0].val,2)} &bull; Bending ${g(c.utils[1].val,2)} &bull; Defl ${g(c.utils[2].val,2)}`
        : sciU
        ? `Shear ${g(c.utils[0].val,2)} &bull; Bending ${g(c.utils[1].val,2)} &bull; LTB ${g(c.utils[2].val,2)} &bull; Defl ${g(c.utils[3].val,2)}`
        : `Shear ${g(c.utils[0].val,2)}   Bending ${g(c.utils[1].val,2)}   ${sec.isBox?'M/Mcx':'LTB'} ${g(c.utils[2].val,2)}   Defl ${g(c.utils[5].val,2)}`}
      ${c.unsupported&&c.unsupported.length?`<br><b>Unsupported exact check(s): ${c.unsupported.length}</b>`:''}</div></div>`;
  // 20 Sep 2026 Beam in Wall (wall.html, js/wall/04-wall-report.js): an optional block between the banner and the brief / report when the page defines reportPrefixHtml(a, c, sec); index.html defines none, so this is '' there
  const prefix=(typeof reportPrefixHtml==='function')? reportPrefixHtml(a,c,sec) : '';

  // 20 Sep 2026 (owner: "I see two briefs; it shall be one brief to avoid
  // confusion"): on the EC3 path the report is the verdict banner followed by
  // the reference software-format brief of js/06-brief-reference.js and nothing
  // else - the load list, the member-forces table, the hover diagrams (loading,
  // shear, moment, deflection, torsion) and every check block live inside the
  // brief in reference software order. The long CED report built below is the
  // BS 5950 output only; its former EC3 branches (SCI restrained /
  // unrestrained blocks, torsion cards, Annex B interaction, EC3 notes, the
  // "Detailed derivation" wrapper) were removed with this change.
  if(S.code==='EC3'){
    let brief;
    if(typeof renderDesignBrief!=='function') brief='<div class="err">Design brief renderer not loaded (js/06-brief-reference.js).</div>';
    else { try{ brief=renderDesignBrief(a,c,sec); } catch(err){ brief=`<div class="err">Design brief could not be rendered: ${err}</div>`; } }
    rep.innerHTML=`
  ${banner}
  ${prefix}
  ${brief}
  <div class="note" style="margin-top:8px;color:#9a8f78">Analysis: 2-node Euler Bernoulli beam elements (direct stiffness); reactions exact, shear/moment by statics, deflection at nodes exact. Section data: SCI P363 section tables. This is a design aid   results to be verified by a competent engineer.</div>`;
    // the section load-line card inside the brief zooms on click / Enter: listeners bound here (no inline handlers in the markup)
    if(typeof sectionViewBindZoom==='function') sectionViewBindZoom(rep);
    return;
  }

  // m-factor formula strings (BS 5950 report)
  const mfRow = sec.isBox
    ? `<div>m<sub>LT</sub></div><div class="formula">SHS   no LTB possible, so m<sub>LT</sub> is not used</div><div class="value">1.000</div><div class="status">cl 4.3.6.1</div>`
    : c.isCant
    ? `<div>m<sub>LT</sub></div><div class="formula">Cantilever ? Table 18</div><div class="value">${g(c.mLT,3)}</div><div class="status">Table 18</div>`
    : `<div>m<sub>LT</sub></div><div class="formula">0.2+(0.15 ${f1(a.Mq,1)}+0.5 ${f1(a.Mh,1)}+0.15 ${f1(a.Mq3,1)})/${f1(Math.abs(a.Mmax),1)}</div><div class="value">${g(c.mLT,3)}</div><div class="status">Table 18</div>`;
  const mxRow = `<div>m<sub>x</sub></div><div class="formula">0.2+(0.1 ${f1(a.Mq,1)}+0.6 ${f1(a.Mh,1)}+0.1 ${f1(a.Mq3,1)})/${f1(Math.abs(a.Mmax),1)} = 0.8 ${f1(a.M24,1)}/${f1(Math.abs(a.Mmax),1)}</div><div class="value">${g(c.mx,3)}</div><div class="status">Table 26</div>`;

  // notes
  const notes=[];
  notes.push('PASS applies only to the implemented member checks and the enabled load combinations. Support bearing, connections and the complete structural system require separate verification.');
  if(S.combos.some(cb=>cb.id!=='c1'&&!cb.sls&&cb.on)) notes.push('Custom ULS factors are used as entered; confirm the complete combination set and favourable/unfavourable actions for the selected standard.');
  if(c.combinationChecks&&c.combinationChecks.length>1) notes.push('BS 5950 m-factors are evaluated separately for every ULS combination. The verdict uses the worst utilisation; the detailed envelope calculation below may have a different governing diagram.');
  if(c.unsupported&&c.unsupported.length) c.unsupported.forEach(n=>notes.push(`<b>NOT COVERED:</b> ${n}`));
  if(c.advisory&&c.advisory.length) c.advisory.forEach(n=>notes.push(`<b>ADVISORY (does not block PASS):</b> ${n}`));
  // BS 5950 notes (the EC3 path returned above with its brief, 20 Sep 2026)
  if(c.shearBuckle) notes.push("d/t &gt; 70e   shear buckling must be checked (cl 4.2.3 / 4.4.5); not covered here (none of the tabulated sections normally reach this limit).");
  if(c.hsNote) notes.push(c.hsNote+".");
  if(c.n>0.1 && !sec.isBox && sec.kind==='I') notes.push(`Significant axial load   the reduced plastic modulus S<sub>rx</sub> uses the standard web-area formula for an equal-flanged I/H section.`);
  if(sec.isBox && S.family==='shs') notes.push(`Closed square section: the BS box-section LTB calculation gives ?<sub>LT</sub> 0 and p<sub>b</sub>=p<sub>y</sub>, so M<sub>b</sub> is governed by M<sub>cx</sub>.`);
  if(sec.isBox && S.family==='rhs') notes.push(`RHS lateral torsional buckling now uses the BS box-section ?<sub>LT</sub>=2.25v(F<sub>b</sub>? <sub>w</sub>) calculation rather than the previous rough Table 15 screen.`);
  if(sec.isBox && S.family==='shs') notes.push(`Strut curve: ${sec.boxType==='CF'?'cold-formed ? curve c (a=5.5)':'hot-finished ? curve a (a=2.0)'} per Table 23, both axes (r<sub>x</sub>=r<sub>y</sub> for a square section). Override in "Robertson const." if a different curve applies.`);
  if(sec.isBox && S.family==='rhs') notes.push(`Strut curve: ${sec.boxType==='CF'?'cold-formed box section ? curve c (a=5.5)':'hot-finished box section ? curve a (a=2.0)'} per Table 23, both axes. Note r<sub>x</sub>?r<sub>y</sub> for a true RHS, so P<sub>c</sub> and P<sub>cy</sub> genuinely differ even though the curve is the same both ways. Override in "Robertson const." if a different curve applies.`);
  if(S.family==='ub') notes.push(`Strut curve (Table 23, rolled I-section): x-x curve ${sec.tf<=40?'a (a=2.0)':'b (a=3.5)'}, y-y curve ${sec.tf<=40?'b (a=3.5)':'c (a=5.5)'} for flange thickness ${sec.tf<=40?'=':'>'}40&nbsp;mm. Override per-axis in "Robertson const." if a different curve applies.`);
  if(S.family==='uc') notes.push(`Strut curve (Table 23, rolled H-section): x-x curve ${sec.tf<=40?'b (a=3.5)':'c (a=5.5)'}, y-y curve ${sec.tf<=40?'c (a=5.5)':'d (a=8.0)'} for flange thickness ${sec.tf<=40?'=':'>'}40&nbsp;mm   one curve lower than a rolled I-section at the same thickness. Override per-axis in "Robertson const." if a different curve applies.`);
  if(sec.kind==='channel') notes.push(`Strut curve: BS&nbsp;5950 directs channel struts to Table 25 (a distinct method from the generic curve a d system). PFC compression is blocked from PASS unless verified Table&nbsp;25/section tables data is implemented.`);
  notes.push("Minor-axis bending M<sub>y</sub> = 0 in this single-plane solver, so m<sub>y</sub> = m<sub>yx</sub> = 1.");

  function classBlockBSfn(){ return sec.isBox ? `
  <div class="section-title">Classification and Properties (BS 5950-1:2000)</div>
  <div class="props">
    <div>Section (${g(sec.mass,1)} kg/m)</div><div>${sname(sec.key)} ${famLabel} [${S.grade}]  D=${g(sec.D,0)} B=${g(sec.B,0)} t=${g(sec.tf,1)} mm</div><div></div><div></div>
    <div>e = v(275/p<sub>y</sub>)</div><div>v(275/${g(a.py,0)})</div><div class="right">${g(c.eps,3)}</div><div></div>
    <div>Flange  b/t = ${g(sec.bT,2)}</div><div>limit ${g(c.cl.flim[0],0)}e=${g(c.cl.flim[0]*c.eps,1)} (Cl.1) / ${g(c.cl.flim[1],0)}e / ${g(c.cl.flim[2],0)}e</div><div class="right">Class ${c.cl.fc}</div><div class="right">${["","Plastic","Compact","Semi-comp","Slender"][c.cl.fc]}</div>
    <div>Web  d/t = ${g(sec.dt,2)}</div><div>limit ${g(c.cl.wlim[0],0)}e=${g(c.cl.wlim[0]*c.eps,1)} (Cl.1) / ${g(c.cl.wlim[1],0)}e / ${g(c.cl.wlim[2],0)}e</div><div class="right">Class ${c.cl.wc}</div><div class="right">${["","Plastic","Compact","Semi-comp","Slender"][c.cl.wc]}</div>
    <div>Class = Fn(b/t,d/t,p<sub>y</sub>,F,M<sub>x</sub>)   Table 12</div><div>${g(sec.bT,2)}, ${g(sec.dt,2)}, ${g(a.py,0)}, ${f1(c.F,1)}, ${f1(a.Mmax,1)}</div><div class="right"></div><div class="right"><b>${c.clsName}</b></div>
  </div>` : `
  <div class="section-title">Classification and Properties (BS 5950-1:2000)</div>
  <div class="props">
    <div>Section (${g(sec.mass,1)} kg/m)</div><div>${sname(sec.key)} ${famLabel} [${S.grade}]</div><div></div><div></div>
    <div>e = v(275/p<sub>y</sub>)</div><div>v(275/${g(a.py,0)})</div><div class="right">${g(c.eps,3)}</div><div></div>
    <div>Flange  b/T = ${g(c.bTBS!=null?c.bTBS:sec.bT,2)}</div><div>b = ${sec.kind==='channel'? 'B (full flange, Figure 5)':'B/2 (rolled I/H, Figure 5)'}; limit 9e = ${g(9*c.eps,1)} (Cl.1) / 10e / 15e</div><div class="right">Class ${c.cl.fc}</div><div class="right">${["","Plastic","Compact","Semi-comp","Slender"][c.cl.fc]}</div>
    ${sec.kind==='channel'
      ? `<div>Web (channel)  d/t = ${g(sec.dt,2)}</div><div>limit 40e = ${g(40*c.eps,1)} (Table 11   flat limit, all classes)</div><div class="right">Class ${c.cl.wc}</div><div class="right">${c.cl.wc<=3?'='+g(40*c.eps,1)+' OK':'Slender'}</div>`
      : `<div>Web (I/H section)  d/t = ${g(sec.dt,2)}</div><div>${c.cl.webCase==='bending+compression'
        ? `bending + compression (Table 11, cl 3.5.5): r<sub>1</sub> = ${g(c.cl.r1,3)}, r<sub>2</sub> = ${g(c.cl.r2,3)}; limits ${g(c.cl.wlim[0],1)}e / ${g(c.cl.wlim[1],1)}e / ${g(c.cl.wlim[2],1)}e`
        : `limit 80e = ${g(80*c.eps,0)} / 100e / 120e`}</div><div class="right">Class ${c.cl.wc}</div><div class="right">${["","Plastic","Compact","Semi-comp","Slender"][c.cl.wc]}</div>`}
    <div>Class = Fn(b/T,d/t,p<sub>y</sub>,F,M<sub>x</sub>,M<sub>y</sub>)</div><div>${g(c.bTBS!=null?c.bTBS:sec.bT,2)}, ${g(sec.dt,2)}, ${g(a.py,0)}, ${f1(c.F,1)}, ${f1(a.Mmax,1)}, 0</div><div class="right"></div><div class="right"><b>${c.clsName}</b></div>
  </div>`; }
  const classBlock = classBlockBSfn();

  function ltbBlockBSfn(){ return sec.isBox ? `
  <div class="section-title smallgap">Lateral Torsional Buckling (Cl. 4.3.6.1${S.family==='rhs'?' / Table 15':''})</div>
  <div class="calc-block">
    <div>L<sub>E</sub> = ${g(c.leK,2)} L</div><div class="formula">${g(c.leK,2)} ${g(S.L,3)} m: ${c.leBasis||''}; end flags ${endsConditionsLine(S)}${c.leMsg? ' &mdash; <b>NOT VERIFIED</b> (no tabulated row: enter L<sub>E</sub>/L)' : ''}</div><div class="value">${g(c.LE/1000,3)} m</div><div class="status">${c.leRow==='user'? 'entered' : c.leMsg? 'blocked' : 'Table 13/14 [verify]'}</div>
    <div>? = L<sub>E</sub>/r<sub>y</sub></div><div class="formula">${f1(c.LE,0)} / ${g(sec.ry*10,1)}</div><div class="value">${f1(c.lam,2)}</div><div></div>
    <div>F<sub>b</sub> = v(S<sub>x</sub> ?'/(AJ))</div><div class="formula">?'=(1-I<sub>y</sub>/I<sub>x</sub>)(1-J/(2.6I<sub>x</sub>))</div><div class="value">${f1(c.phiB,3)}</div><div></div>
    <div> <sub>w</sub></div><div class="formula">${c.cl.cls<=2?'Class 1/2 ? 1.0':'Z<sub>x</sub>/S<sub>x</sub>'}</div><div class="value">${g(c.betaW,3)}</div><div></div>
    <div>?<sub>LT</sub> = 2.25v(F<sub>b</sub>? <sub>w</sub>)</div><div class="formula">2.25v(${f1(c.phiB,3)} ${f1(c.lam,2)} ${g(c.betaW,3)})</div><div class="value">${f1(c.lamLT,2)}</div><div></div>
    <div>p<sub>b</sub> = Fn(p<sub>y</sub>, ?<sub>LT</sub>)</div><div class="formula">?<sub>L0</sub> = ${f1(c.lamL0,1)}</div><div class="value">${f1(c.pb,1)} N/mm </div><div class="status">Annex B</div>
    <div>M<sub>b</sub> = ${c.cl.cls<=2?'S':'Z'}<sub>x</sub> p<sub>b</sub> = M<sub>cx</sub></div><div class="formula">${g(c.cl.cls<=2?sec.Sx:sec.Zx,1)} ${f1(c.pb,1)}</div><div class="value">${f1(c.Mb,2)} kN m</div>${st(c.ltbUtil<=1,'OK')}
  </div>` : `
  <div class="section-title smallgap">Lateral Torsional Buckling   M<sub>b</sub> (Cl. 4.3 / Annex B)</div>
  <div class="calc-block">
    <div>L<sub>E</sub> = ${g(c.leK,2)} L</div><div class="formula">${g(c.leK,2)} ${g(S.L,3)} m: ${c.leBasis||''}; end flags ${endsConditionsLine(S)}${c.leMsg? ' &mdash; <b>NOT VERIFIED</b> (no tabulated row: enter L<sub>E</sub>/L)' : ''}; strut lengths L<sub>cr,x</sub> = ${g(c.LcrX/1000,3)} m, L<sub>cr,y</sub> = ${g(c.LcrY/1000,3)} m (${c.lcrBasis||''})</div><div class="value">${g(c.LE/1000,3)} m</div><div class="status">${c.leRow==='user'? 'entered' : c.leMsg? 'blocked' : 'Table 13/14 [verify]'}</div>
    <div>? = L<sub>E</sub> / r<sub>y</sub></div><div class="formula">${f1(c.LE,0)} / ${g(sec.ry*10,1)}</div><div class="value">${f1(c.lam,2)}</div><div></div>
    <div>v = 1/[1+0.05(?/x) ]<sup> </sup></div><div class="formula">x = ${g(sec.x,1)} (torsional index)</div><div class="value">${g(c.v,3)}</div><div class="status">N=0.5</div>
    <div> <sub>w</sub></div><div class="formula">${c.cl.cls<=2?'Class 1/2 ? 1.0':'Z<sub>x</sub>/S<sub>x</sub>'}</div><div class="value">${g(c.betaW,3)}</div><div></div>
    <div>?<sub>LT</sub> = u v ? v <sub>w</sub></div><div class="formula">${g(sec.u,3)} ${g(c.v,3)} ${f1(c.lam,2)} v${g(c.betaW,2)}</div><div class="value">${f1(c.lamLT,2)}</div><div></div>
    <div>p<sub>b</sub> = Fn(p<sub>y</sub>, ?<sub>LT</sub>)</div><div class="formula">?<sub>L0</sub> = ${f1(c.lamL0,1)} ${c.lamLT<=c.lamL0?'(?<sub>LT</sub>=?<sub>L0</sub> ? p<sub>b</sub>=p<sub>y</sub>)':'(Perry, Annex B.2)'}</div><div class="value">${f1(c.pb,1)} N/mm </div><div class="status">Annex B</div>
    <div>M<sub>b</sub> = ${c.cl.cls<=2?'S':'Z'}<sub>x</sub> p<sub>b</sub> = M<sub>cx</sub></div><div class="formula">${g(c.cl.cls<=2?sec.Sx:sec.Zx,1)} ${f1(c.pb,1)}</div><div class="value">${f1(c.Mb,2)} kN m</div><div></div>
    <div>m<sub>LT</sub>M<sub>x</sub> = M<sub>b</sub> (cl 4.3.6.2)</div><div class="formula">${g(c.mLT,3)} ${f1(c.Mx,1)} / ${f1(c.Mb,2)}</div><div class="value">${g(c.ltbUtil,3)}</div>${st(c.ltbUtil<=1,'OK')}
  </div>`; }
  const ltbBlock = ltbBlockBSfn();

  const pvFormula = sec.isBox? (sec.D===sec.B? 'P<sub>v</sub>=0.6 p<sub>y</sub> A D/(D+B)=0.6 p<sub>y</sub> A/2' : 'P<sub>v</sub>=0.6 p<sub>y</sub> A D/(D+B)') : 'P<sub>v</sub>=0.6 p<sub>y</sub> t D';

  rep.innerHTML = `
  ${banner}
  ${prefix}
  <div class="report-head">
    <div>
      <h2>Member Loading and Member Forces</h2>
      <div class="meta">${sname(sec.key)} ${famLabel} &nbsp; &nbsp; ${gradeTxt} &nbsp; &nbsp; L = ${g(S.L)} m</div>
      ${endsLine}
      <div class="loadlist">${loadLines}</div>
    </div>
    ${sectionView}
  </div>

  <div class="diagcard">
    <div class="dt">Loading</div>
    ${beamDiagram(a)}
    <div class="note" style="margin-left:0">Reactions (governing-moment combo, ${a.governM.combo.label}; R upward positive, M sagging positive): ${reactLine}</div>
    ${upliftLines}
  </div>

  ${combosBlock}

  <table class="force-table">
    <thead>
      <tr><th class="table-title" colspan="4">Member Force Envelope   ULS combinations&nbsp;&nbsp;|&nbsp;&nbsp;Deflection   SLS combinations</th></tr>
      <tr><th>Axial Force (kN)</th><th>Max Shear F<sub>v</sub> (kN)</th><th>Max Moment M<sub>x</sub> (kN m @ m)</th><th>Max Deflection (mm @ m)</th></tr>
    </thead>
    <tbody><tr>
      <td class="num">${f1(c.F,3)} ${c.F>=0?'C':'T'}</td>
      <td class="num">${f1(a.Vmax,2)}</td>
      <td class="num">${f1(a.Mmax,2)} @ ${g(a.Mpos,2)}</td>
      <td class="num">${f1(a.dmax,2)} @ ${g(a.dpos,2)}</td>
    </tr></tbody>
  </table>

  <div class="diagcard">
    <div class="diagrow">
      <div><div class="dt">Shear force (kN)</div>${plot(a.diag.xs,a.diag.V,{color:'#1d4ed8',fill:'#bcd0f7',unit:'kN',fmt:v=>f1(v,2)})}</div>
      <div><div class="dt">Bending moment (kN m)</div>${plot(a.diag.xs,a.diag.M,{color:'#b91c1c',fill:'#f3c2c2',unit:'kN m',flip:true,fmt:v=>f1(v,2)})}</div>
    </div>
    <div class="diagrow" style="margin-top:6px">
      <div><div class="dt">Deflection (mm, SLS)</div>${plot(a.diag.dx,a.diag.dw,{color:'#166534',fill:'#bfe3cb',unit:'mm',fmt:v=>f1(v,2)})}</div>
      <div style="display:flex;align-items:center;padding:8px 4px"><div class="note" style="margin:0">Bending-moment diagram plotted on the tension (sagging-down) side. Shear and deflection plotted to true sign (down = below the axis).</div></div>
    </div>
  </div>

  <div class="calcs-start"></div>
  ${classBlock}

  <div class="section-title smallgap">Local Capacity Check (Cl. 4.2)</div>
  <div class="calc-block">
    <div>F<sub>vx</sub> / P<sub>vx</sub></div><div class="formula">${f1(c.Fv,2)} / ${f1(c.Pv,1)} &nbsp;(${pvFormula})</div><div class="value">${g(c.Fv/c.Pv,3)}</div>${st(c.lowShear,'Low Shear','High Shear')}
    <div>M<sub>cx</sub> = p<sub>y</sub> S<sub>x</sub> = 1.2 p<sub>y</sub> Z<sub>x</sub></div><div class="formula">${g(a.py,0)} ${g(sec.Sx,1)} = 1.2 ${g(a.py,0)} ${g(sec.Zx,1)}</div><div class="value">${f1(c.Mcx,2)} kN m</div><div></div>
    <div>A<sub>e</sub> = K<sub>e</sub> A<sub>net</sub> = A<sub>g</sub></div><div class="formula">${g(c.Ke,2)} ${g(c.Anet/1e2,1)} = ${g(c.Ag/1e2,1)}</div><div class="value">${g(c.Ae/1e2,1)} cm </div><div></div>
    <div>P<sub>z</sub> = A<sub>e</sub> p<sub>y</sub></div><div class="formula">${g(c.Ae/1e2,1)} ${g(a.py,0)}</div><div class="value">${f1(c.Pz,1)} kN</div><div></div>
    <div>n = F / P<sub>z</sub></div><div class="formula">${f1(c.F,2)} / ${f1(c.Pz,1)}</div><div class="value">${g(c.n,3)}</div>${st(c.n<=1,'OK')}
    <div>S<sub>rx</sub> = Fn(S<sub>x</sub>, n)</div><div class="formula">${g(sec.Sx,1)}${c.n>0.02?` reduced for n`:`, n 0`}</div><div class="value">${g(c.Srx/1e3,1)} cm </div><div></div>
    <div>M<sub>rx</sub> = min(M<sub>cx</sub>, axial-reduced resistance)</div><div class="formula">including the applicable shear reduction and elastic cap</div><div class="value">${f1(c.Mrx,2)} kN m</div><div></div>
    <div>(M<sub>x</sub>/M<sub>rx</sub>) + (M<sub>y</sub>/M<sub>ry</sub>)</div><div class="formula">(${f1(c.Mx,2)}/${f1(c.Mrx,2)}) + 0</div><div class="value">${g(c.localUtil,3)}</div>${st(c.localUtil<=1,'OK')}
  </div>

  <div class="section-title smallgap">Equivalent Uniform Moment Factors</div>
  <div class="calc-block">
    ${mfRow}
    ${mxRow}
    <div>m<sub>y</sub>, m<sub>yx</sub></div><div class="formula">M<sub>y</sub> = 0 (single-plane bending)</div><div class="value">1.000</div><div class="status">Table 26</div>
  </div>

  ${ltbBlock}

  <div class="section-title smallgap">Simplified Buckling Approach (Cl. 4.8.3.3.1)</div>
  <div class="calc-block">
    <div>p<sub>y</sub> Z<sub>x</sub></div><div class="formula">${g(a.py,0)} ${g(sec.Zx,1)}</div><div class="value">${f1(c.pyZx,2)} kN m</div><div></div>
    <div>F/P<sub>c</sub> + m<sub>x</sub> M<sub>x</sub>/(p<sub>y</sub>Z<sub>x</sub>)</div><div class="formula">${f1(c.Fc,1)}/${f1(c.Pc,0)} + ${g(c.mx,3)} ${f1(c.Mx,2)}/${f1(c.pyZx,2)}</div><div class="value">${g(c.u1,3)}</div>${st(c.u1<=1,'OK')}
    <div>F/P<sub>cy</sub> + m<sub>LT</sub> M<sub>LT</sub>/M<sub>b</sub></div><div class="formula">${f1(c.Fc,1)}/${f1(c.Pcy,0)} + ${g(c.mLT,3)} ${f1(c.Mx,2)}/${f1(c.Mb,2)}</div><div class="value">${g(c.u2,3)}</div>${st(c.u2<=1,'OK')}
  </div>

  <div class="section-title smallgap">Deflection Check (SLS &mdash; ${a.governD.combo.label})</div>
  <div class="calc-block">
    <div>${c.deflCant? 'Tip deflection (vertically free end, relative to the held end)' : 'w (governing span utilisation)'}</div><div class="formula">@ x = ${g(a.deflection?a.deflection.dpos/1000:a.dpos,2)} m</div><div class="value">${f1(c.dmax,1)} mm</div><div></div>
    <div>Limit = ${c.deflCant? 'L/' : 'span/'}${g(c.divisor,0)}${a.deflection&&a.deflection.abs!=null? ', capped at '+f1(a.deflection.abs,1)+' mm (absolute)' : ''}</div><div class="formula">${g(c.span,0)}/${g(c.divisor,0)} = ${f1(a.deflection?a.deflection.limSpan:c.dlimit,1)} mm${c.deflAbsGoverns? ' &gt; absolute limit '+f1(a.deflection.abs,1)+' mm, which governs' : ''}${c.deflCant? ' (vertically free end: L/'+g(c.divisor,0)+' per UK NA to EN 1993-1-1 Table NA.2 [verify], cantilever row)' : ''}</div><div class="value">${f1(c.dmax,1)} ${c.defOk?'&lt;':'&gt;'} ${f1(c.dlimit,1)} mm</div>${st(c.defOk,'OK')}
  </div>

  <div class="note">${notes.map(n=>'  '+n).join('<br>')}</div>
  <div class="note" style="margin-top:8px;color:#9a8f78">Analysis: 2-node Euler Bernoulli beam elements (direct stiffness); reactions exact, shear/moment by statics, deflection at nodes exact. Section data: SCI P363 section tables. This is a design aid   results to be verified by a competent engineer.</div>`;
  // the section load-line card zooms on click / Enter: listeners bound here
  // (no inline handlers in the markup)
  if(typeof sectionViewBindZoom==='function') sectionViewBindZoom(rep);
}


/* ==== beam-v03 module 15 ==== */
/* ===========================================================================
   6b. REFERENCE SOFTWARE-FORMAT DESIGN BRIEF  (EC3 path)
   ---------------------------------------------------------------------------
   renderDesignBrief(a, c, sec) -> HTML string.
   Pure: reads the analysis result a (analyse()), the check result c
   (checks(a)), the normalised section sec (= a.sec) and the app state S;
   touches no DOM. Block order, wording and the three-column layout follow
   docs/BRIEF_MAPPING.md, i.e. a reference software 2025 "MasterSteel Beam Design
   to EN 1993-1-1" printout. Every number printed is a field of a / c / sec
   or is derived by one of the msb* helpers below, each of which only
   re-expresses a value the check engine already used (ratio, unit change,
   back-substitution); no resistance is recomputed here (V_pl.z.Rd and
   M_c.z.Rd come from the engine's c.ax block). Where the engine has no
   value the tag column prints "not evaluated".
   Blocking messages (c.unsupported) appear as red "NOT VERIFIED" rows in
   the block of the check they concern (a message no block claims is printed
   before the deflection block) and again in the verdict footer.
   20 Sep 2026 (owner: one brief only, reference software order, every block on its
   trigger): this brief IS the EC3 report (js/06-render.js prints the verdict
   banner and then this string alone). Block sequence, audited against the
   reference software printouts of docs/owner-cases/reference/COMPARISON_LOG.md:
   Title (+ "Includes Design for Torsion ..." line when torsion is active) ->
   Member Loading and Member Forces (load list, forces table, loading sketch
   and hover diagrams V / M / delta / T, combination table when > 1 case) ->
   Classification and Effective Area -> Shear Capacity Check -> Local Capacity
   Check (or Moment Capacity Check M.c.y.Rd) -> Web Transverse Forces (beam-v03
   addition, when the engine ran it) -> Compression Resistance N.b.Rd (N_Ed > 0)
   -> Equivalent Uniform Moment Factor(s) -> Lateral Buckling Check M.b.Rd
   (+ Lateral Restraint Portions with intermediate restraints) -> Buckling
   Resistance (Axial with Moments brief) -> Torsion Design / Torsion Bending
   Design @ x / Torsion Shear Design @ x (torsion active) -> Deflection Check
   (+ the "Torq in Case n" twist line with torsion) -> unity bar -> verdict.
   A reference software line beam-v03 has no value for prints its label with
   "n/a - not evaluated by beam-v03" once; nothing is recomputed here.
   20 Sep 2026 torsion + N/Mz (engine: combined torsion with N_Ed and / or an
   imposed M_z verified as Eurocode advises, the former "not implemented as one
   interaction" block gone): the Local Capacity / Moment Capacity block prints
   the EN 1993-1-1 6.2.7(5) elastic yield criterion (6.1) row with its stress
   build-up and the section points (c.tor.elastic, every torsion case); the
   Torsion Design block prints M_z,tot = M_z,Ed + phi.M_y, the (A.1) rows with
   M_z,tot and C_mz = 1.0 with its reason, the basis text (tor.combinedBasis)
   and - N_Ed > 0 - the superposition of Eq 6.62 and (A.1) as an information-only
   row (tor.superposition, class ms-advrow); the Buckling Resistance U_M.z row
   and the C_mz row print M_z,Ed as imposed + twist = total; the unity bar
   carries the "Yield 6.1" cell.
   =========================================================================== */

/* ---- formatting (docs/BRIEF_MAPPING.md section 3) ---- */
function msbR(v){ return f1(v,3); }                       // ratio / chi / k / lambda-bar
function msbKNm(v){ return f1(v,3); }                     // kN.m
function msbKN(v){ return f1(v,3); }                      // kN
function msbMM(v){ return f1(v,2); }                      // mm
function msbM(v){ return g(v,3); }                        // m, trimmed
function msbInt(v){ return g(v,0); }
function msbDash(v,fmt){ return (v==null||!isFinite(v))? '&mdash;' : fmt(v); }
function msbEsc(s){ return String(s==null?'':s); }
function msbSecName(key){ return String(key||'').replace(/\s*x\s*/g,' x '); }   // "457 x 191 x 89", "150 x 150 x 6.3"

/* ---- pure DERIVE helpers (docs/BRIEF_MAPPING.md section 8) ---- */
// shear at a station of one combination's own diagram, kN
function msbEndShear(fb,x){ return interpAt(fb.xs,fb.V,x)/1000; }
// bending moment at a station, kN.m (signed, sagging positive)
function msbMomentAt(fb,x){ return interpAt(fb.xs,fb.M,x)/1e6; }
// reference software C1 inputs of a portion: M1, M2 (|M2| >= |M1|), Mo (mid-portion
// moment above the chord), psi = M1/M2, mu = Mo/M2 capped at 300; plus the
// largest |M| on the grid inside the portion and its position
function msbPortionMoments(fb,xa,xb){
  const Ma=msbMomentAt(fb,xa+1e-4), Mb=msbMomentAt(fb,xb-1e-4);
  const [M1,M2]= Math.abs(Mb)>=Math.abs(Ma)? [Ma,Mb] : [Mb,Ma];
  const Mmid=mAtStation(fb,(xa+xb)/2,xa,xb)/1e6;   // larger side of a jump (in-span couple), as c1Inputs
  const Mo=Mmid-(M1+M2)/2;
  const psi= Math.abs(M2)>1e-9? M1/M2 : 0;
  const mu= Math.abs(M2)>1e-9? Math.min(Mo/M2,300) : 300;
  let Mmax=0, xmax=xa;
  fb.xs.forEach((x,i)=>{ if(x>=xa-1e-6&&x<=xb+1e-6&&Math.abs(fb.M[i])>Math.abs(Mmax)){ Mmax=fb.M[i]; xmax=x; } });
  return {M1,M2,Mo,psi,mu,Mmax:Mmax/1e6,xmax};
}
// Table B.3 parameters exactly as cmTableB3() reads them from the analysis
function msbCmB3Params(a){
  const M0=a.M0end, ML=a.MLend, Ms=a.Mh;
  const Mh= Math.abs(M0)>=Math.abs(ML)? M0 : ML;
  const Mo= Math.abs(M0)>=Math.abs(ML)? ML : M0;
  const psi= Math.abs(Mh)>1e-9? Math.max(-1,Math.min(1,Mo/Mh)) : 1;
  const useAlphaS = Math.abs(Ms)<=Math.abs(Mh)+1e-12;
  const alphaS = useAlphaS? (Math.abs(Mh)>1e-9? Ms/Mh : 0) : null;
  const alphaH = useAlphaS? null : Mh/Ms;
  return {Mh,Ms,psi,alphaS,alphaH};
}
// Table B.3 expression text from the annexB2 cmLabel prefix
function msbCmB3Form(label,alphaS){
  const s=msbEsc(label);
  const hasH=/alpha;<sub>h<\/sub>/.test(s), hasS=/alpha;<sub>s<\/sub>/.test(s);
  if(s.indexOf('linear end-moment diagram')===0) return 'Max(0.6+0.4&psi;, 0.4)';
  if(s.indexOf('negligible moment')===0) return '1.0';
  if(s.indexOf('C_m = 1')===0) return '1.0 (diagram outside Table B.3)';
  if(s.indexOf('mixed loading')===0) return 'max(uniform, concentrated)';
  if(s.indexOf('uniform load diagram')===0){
    if(hasH) return '0.95+0.05&alpha;<sub>h</sub>';
    if(hasS) return (alphaS!=null&&alphaS<0)? 'Table B.3 (&alpha;<sub>s</sub> &lt; 0 row)' : '0.2+0.8&alpha;<sub>s</sub>';
    return '0.95 (M<sub>h</sub> = 0)';
  }
  if(s.indexOf('concentrated load diagram')===0){
    if(hasH) return '0.90+0.10&alpha;<sub>h</sub>';
    if(hasS) return (alphaS!=null&&alphaS<0)? 'Table B.3 (&alpha;<sub>s</sub> &lt; 0 row)' : '0.2+0.8&alpha;<sub>s</sub>';
    return '0.90 (M<sub>h</sub> = 0)';
  }
  return 'Table B.3';
}
// reference software C1 tag from the closed-form label / route
function msbC1Tag(label,route){
  const s=msbEsc(label);
  if(route==='override'||/user override/.test(s)) return 'User';
  if(route==='sn006a'||route==='cantilever'||/cantilever/.test(s)) return 'Cantilever';
  if(/uniformly distributed/.test(s)) return 'Uniform';
  if(/central point load/.test(s)) return 'Point';
  if(/linear end-moment gradient/.test(s)) return 'Not Loaded';
  if(/Serna/.test(s)) return 'Serna';
  if(/negligible/.test(s)) return 'Uniform';
  return '';
}
// Euler load about one axis, kN: pi^2 E I / L^2 (I cm4, L mm)
function msbNcr(E,I_cm4,L_mm){ return Math.PI*Math.PI*E*I_cm4*1e4/(L_mm*L_mm)/1000; }
// characteristic resistances used by the cl 6.3.3 ratios
function msbNRk(sec,fy){ return sec.A*100*fy/1000; }           // kN
function msbMyRk(Wy,fy){ return Wy*fy/1e6; }                   // kN.m (Wy mm3)
function msbKc(C1){ return (typeof kcFromC1==='function')? kcFromC1(C1).kc : Math.max(Math.min(1/Math.sqrt(Math.max(C1,1e-6)),1),1/Math.sqrt(2.76)); }   // NA 2.18 with the Table 6.6 floor 0.60
// cl 6.2.8(3) shear reduction factor
function msbRhoShear(V,Vpl){ return Math.min(Math.pow(2*V/Math.max(Vpl,1e-9)-1,2),1); }
// reduced plastic modulus back-substituted from M_N,Rd, cm3
function msbWplN(MN_kNm,fy){ return MN_kNm*1e3/fy; }
// Phi_LT of cl 6.3.2.3 (the expression the check engine evaluates)
function msbPhiLT(lam,alphaLT){ return 0.5*(1+alphaLT*(lam-0.4)+0.75*lam*lam); }
// 1-based index of a combination among the analysed ULS (or SLS) combinations
// (the list of analyse() when `a` is given, else the enabled entries of S.combos)
function msbCaseIndex(combo,sls,a){
  const list= a? (sls? a.slsResults : a.ulsResults).map(r=>r.combo) : (S.combos||[]).filter(cb=>cb.on && (sls? cb.sls : !cb.sls));
  let i=list.indexOf(combo);
  if(i<0 && combo) i=list.findIndex(cb=>cb.label===combo.label);
  return i<0? null : i+1;
}
function msbCaseLabel(combo,sls,a){
  const n=msbCaseIndex(combo,sls,a);
  return (n!=null? String(n) : '?')+(combo&&combo.label? ' ('+combo.label+')' : '');
}
// [1,2,4,5,6] -> "1-2, 4-6"
function msbCaseRanges(list){
  const out=[]; let i=0;
  while(i<list.length){ let j=i; while(j+1<list.length && list[j+1]===list[j]+1) j++; out.push(list[i]===list[j]? String(list[i]) : list[i]+'-'+list[j]); i=j+1; }
  return out.join(', ');
}
function msbAxialTag(sec,eps){ return ((S.axial||0)>0 && sec.dt>42*eps)? '(Axial: Slender web)' : '(Axial: Non-Slender)'; }
function msbMaxExclDeflection(utils){
  const v=(utils||[]).filter(u=>u.name!=='Deflection').map(u=>+u.val).filter(x=>isFinite(x));
  return v.length? Math.max(...v) : 0;
}
// governing portion [xa, xb] in mm (docs/BRIEF_MAPPING.md section 7.1)
function msbPortion(a,LT){
  const L=a.L;
  if(LT && LT.spanGoverns && LT.spanGov) return [LT.spanGov.a, LT.spanGov.b];
  const pts=(LT && LT.vPoints)||[];
  if(pts.length>=3 && LT.modePeakX!=null){
    for(let i=0;i<pts.length-1;i++) if(LT.modePeakX>=pts[i]-1e-6 && LT.modePeakX<=pts[i+1]+1e-6) return [pts[i],pts[i+1]];
  }
  return [0,L];
}
// which brief block a blocking message belongs to
function msbBlockFor(msg){
  const s=msbEsc(msg);
  if(/^Hold-down/.test(s)) return 'forces';
  if(/Web transverse forces|bearing stiffener|EN 1993-1-5 clause 6|EN 1993-1-5 9\.4/i.test(s)) return 'web';
  if(/torsion|Torsion|Eccentric loads|Annex A|k_alpha|Torsional constants|eccentric load/i.test(s)) return 'torsion';
  if(/Class 4|effective-section|Effective-area|slender/i.test(s)) return 'class';
  if(/shear/i.test(s)) return 'local';
  if(/PFC under axial|torsional-flexural|6\.3\.1\.4/.test(s)) return 'compression';
  if(/reduced moment resistance under axial|channel sections/i.test(s)) return 'local';
  if(/critical moment|LTB|eigensolver|mesh|Cantilever|lateral/i.test(s)) return 'ltb';
  return 'general';
}
/* ---- [20 Sep 2026 torsion + N/Mz] elastic yield criterion (6.1) rows, c.tor.elastic ----
   EN 1993-1-1 cl 6.2.7(5) with 6.2.1(5) Eq (6.1), (sigma_x/f_yd)^2 + 3(tau/f_yd)^2 <= 1,
   evaluated by torsionElasticOpen / torsionElasticBox (js/checks/eurocode-checks.js)
   at every station of every ULS combination and at the section points: open
   sections P1 flange tip, P2 / P3 the web-flange junction (flange / web side), P4
   the web mid-depth; hollow sections (6.2.7(7), warping neglected) the corner,
   the web mid-depth and the flange mid-width. msbElasticTerms lists, point by
   point, the formula term of every stress component the engine stored (a term
   the engine sets to zero for the section - W_n2 / S_w2 / S_w3 of a doubly
   symmetric I/H - is not printed); the printed sum is the engine's sigmaX / tau
   and the value its u. Nothing is recomputed here.
   20 Sep 2026 review: P2 carries the flange shear flow V S_f/(2 I_y t_f) (I/H, half flange) or V S_f/(I_y t_f)
   (channel, whole flange) and - channel - E S_w2 phi'''/t_f at the junction; the channel point P1b (flange point
   where W_n = 0, S_w1) prints its own terms; the box corner carries V Q_c/(I_y t) and the box web mid-depth
   V Q_m/(I_y t) (the closed-section shear flow, Eq 6.20). */
function msbElasticTerms(p,el){
  const chan=!!(el.geom && el.geom.chan), box=!!el.box, mz= el.MzImp>1e-9? 'M<sub>z,tot</sub>' : 'M<sub>z</sub>';
  const sN={l:'N/A',k:'sigmaN'}, sMyW={l:'M<sub>y</sub>/W<sub>el.y</sub>',k:'sigmaMy'}, sMzW={l:mz+'/W<sub>el.z</sub>',k:'sigmaMz'}, sMzWeb={l:mz+'y<sub>w</sub>/I<sub>z</sub>',k:'sigmaMz'};
  const sW0={l:'EW<sub>n0</sub>&phi;&Prime;',k:'sigmaW'}, sW2={l:'EW<sub>n2</sub>&phi;&Prime;',k:'sigmaW'};
  const tTf={l:'Gt<sub>f</sub>&phi;&prime;',k:'tauT'}, tTw={l:'Gt<sub>w</sub>&phi;&prime;',k:'tauT'}, tT={l:'T<sub>Ed</sub>/W<sub>t</sub>',k:'tauT'};
  const sMyFm={l:'M<sub>y</sub>(D/2&minus;t<sub>f</sub>/2)/I<sub>y</sub>',k:'sigmaMy'};
  if(box){
    if(/^corner/.test(p.point)) return {sig:[sN,sMyW,sMzW], tau:[tT,{l:'VQ<sub>c</sub>/(I<sub>y</sub>t)',k:'tauV'}]};
    if(/^web mid/.test(p.point)) return {sig:[sN,sMzW], tau:[tT,{l:'VQ<sub>m</sub>/(I<sub>y</sub>t)',k:'tauV'}]};
    return {sig:[sN,sMyW], tau:[tT]};
  }
  if(/^P1b/.test(p.point)) return {sig:[sN,sMyFm,{l:mz+'y<sub>1</sub>/I<sub>z</sub>',k:'sigmaMz'}], tau:[{l:'Vs<sub>1</sub>t<sub>f</sub>(D&minus;t<sub>f</sub>)/(2I<sub>y</sub>t<sub>f</sub>)',k:'tauV'},tTf,{l:'ES<sub>w1</sub>&phi;&#8244;/t<sub>f</sub>',k:'tauW'}]};
  if(/^P1/.test(p.point)) return {sig:[sN,sMyW,sMzW,sW0], tau:[tTf]};
  if(/^P2/.test(p.point)) return {sig:[sN,sMyFm,sMzWeb].concat(chan? [sW2] : []),
                                  tau:[{l:chan? 'VS<sub>f</sub>/(I<sub>y</sub>t<sub>f</sub>)' : 'VS<sub>f</sub>/(2I<sub>y</sub>t<sub>f</sub>)',k:'tauV'},tTf,{l:chan? 'ES<sub>w2</sub>&phi;&#8244;/t<sub>f</sub>' : 'ES<sub>w1</sub>&phi;&#8244;/t<sub>f</sub>',k:'tauW'}]};
  if(/^P3/.test(p.point)) return {sig:[sN,{l:'M<sub>y</sub>(D/2&minus;t<sub>f</sub>)/I<sub>y</sub>',k:'sigmaMy'},sMzWeb].concat(chan? [sW2] : []),
                                  tau:[{l:'VS<sub>f</sub>/(I<sub>y</sub>t<sub>w</sub>)',k:'tauV'},tTw].concat(chan? [{l:'ES<sub>w2</sub>&phi;&#8244;/t<sub>w</sub>',k:'tauW'}] : [])};
  return {sig:[sN,sMzWeb], tau:[{l:'VS<sub>max</sub>/(I<sub>y</sub>t<sub>w</sub>)',k:'tauV'},tTw].concat(chan? [{l:'ES<sub>w3</sub>&phi;&#8244;/t<sub>w</sub>',k:'tauW'}] : [])};
}
// "N/A + M_y/W_el.y + ... = a + b + ... = total" (one term: "label = value"); total = the engine's sigmaX / tau
function msbElasticSum(p,terms,total){
  if(terms.length===1) return terms[0].l+' = '+f1(total,2);
  return terms.map(t=>t.l).join(' + ')+' = '+terms.map(t=>f1(p[t.k],2)).join(' + ')+' = '+f1(total,2);
}
// the two rows of the Local Capacity block: the governing station / point with its stress build-up reference software-style,
// then the compact values of every section point at that station (class ms-pts)
function msbElasticRows(T){
  const el=T.elastic, pts=el.points||[], gov=pts.find(p=>p.point===el.point)||el, fyd=el.fy;
  const imp=el.MzImp>1e-9;
  const stn='@ x = '+msbM(el.x/1000)+' m ('+msbEsc(el.combo)+'): M<sub>y</sub> = '+msbKNm(el.My)+' kN.m, V = '+msbKN(el.V)+' kN, '+
    (el.box? 'T = '+msbKNm(el.T||0)+' kN.m, &phi; = '+f1(el.phi,5)+' rad'
           : '&phi; = '+f1(el.phi,5)+' rad, &phi;&prime; = '+Math.abs(el.p1).toExponential(3)+' rad/mm, &phi;&Prime; = '+Math.abs(el.p2).toExponential(3)+' rad/mm&sup2;, &phi;&#8244; = '+Math.abs(el.p3).toExponential(3)+' rad/mm&sup3;')+
    '; '+(imp? 'M<sub>z,tot</sub> = M<sub>z.Ed</sub> + &phi;M<sub>y</sub> = '+msbKNm(el.MzImp)+' + '+msbKNm(el.MzTwist)+' = '+msbKNm(el.MzTot) : 'M<sub>z</sub> = &phi;M<sub>y</sub> = '+msbKNm(el.MzTwist))+' kN.m';
  const tm=msbElasticTerms(gov,el);
  const vals=stn+'; '+msbEsc(gov.point)+': &sigma;<sub>x</sub> = '+msbElasticSum(gov,tm.sig,gov.sigmaX)+' N/mm&sup2;; &tau; = '+msbElasticSum(gov,tm.tau,gov.tau)+' N/mm&sup2;; '+
    '('+f1(gov.sigmaX,2)+'/'+msbInt(fyd)+')&sup2; + 3('+f1(gov.tau,2)+'/'+msbInt(fyd)+')&sup2; =';
  // 20 Sep 2026 review: the tag says whether (6.1) binds the verdict (elasticBindingPolicy: Class 3, or an open Class 1/2
  // section with N_Ed) or is printed for information (Class 1/2 with the plastic route of 6.2.7(6) / 6.2.7(7)-(9))
  const bind=!!el.binding;
  const tag= bind? msbWarn(el.u<=1.0001)+' 6.2.7(5)' : (el.u<=1.0001? '&le; 1' : '<span class="ms-warn">&gt; 1</span>')+' information 6.2.7(5)';
  let h=msbRow('Elastic yield criterion with torsion (cl 6.2.7(5), Eq 6.1) @ '+msbM(el.x/1000)+' m, '+msbEsc(el.point), vals, msbR(el.u), tag);
  const cmp=(p)=>{ const t=msbElasticTerms(p,el); return msbEsc(p.point)+': &sigma;<sub>x</sub> = '+(t.sig.length>1? t.sig.map(q=>f1(p[q.k],2)).join(' + ')+' = ' : '')+f1(p.sigmaX,2)+', &tau; = '+(t.tau.length>1? t.tau.map(q=>f1(p[q.k],2)).join(' + ')+' = ' : '')+f1(p.tau,2)+' &rarr; '+msbR(p.u)+(p.point===el.point? ' (governs)' : ''); };
  const gm=el.geom||{};
  h+=msbRow('&sigma;<sub>x</sub>, &tau;, (6.1) at the '+pts.length+' section points', pts.map(cmp).join('; ')+(el.box? ' (Q<sub>c</sub> = '+g((gm.Qc||0)/1e3,1)+', Q<sub>m</sub> = '+g((gm.Qm||0)/1e3,1)+' cm&sup3;: the V<sub>Ed</sub> shear flow round the closed mid-line, Eq 6.20)' : gm.chan? ' (y<sub>w</sub> = web back from the minor axis = '+msbMM(gm.yWeb)+' mm; P1b at s<sub>1</sub> = '+msbMM(gm.s1)+' mm from the toe, y<sub>1</sub> = '+msbMM(gm.y1b)+' mm; S<sub>w1</sub>, S<sub>w2</sub>, S<sub>w3</sub> = P385 Table A.2 points 1, 2, 3)' : ' (y<sub>w</sub> = t<sub>w</sub>/2 = '+msbMM(gm.yWeb)+' mm; S<sub>f</sub>/2 = half flange one side of the web)'),
    '', 'f<sub>y</sub>/&gamma;<sub>M0</sub> = '+msbInt(fyd), 'ms-pts');
  if(el.bindingBasis) h+=msbRow('Elastic verification (6.1): '+(bind? 'verdict-binding' : 'information only'), msbEsc(el.bindingBasis), '', bind? '6.2.7(5)' : '6.2.7(5)-(6)', 'ms-basis');
  return h;
}
// [20 Sep 2026 torsion + N/Mz] the grid row of the largest |phi.M_y| (the M_z,tot line of the Torsion block), a lookup in c.tor.grids
function msbMzTwistRow(T){
  let best=null;
  (T.grids||[]).forEach(g2=>g2.rows.forEach(r2=>{ if(!best || r2.Mz>best.Mz) best=Object.assign({combo:g2.combo.label},r2); }));
  return best;
}

/* ---- row / block builders ---- */
function msbWarn(ok){ return ok? 'OK' : '<span class="ms-warn">Warning</span>'; }
function msbRow(label,vals,res,tag,cls){
  return '<div class="ms-row'+(cls? ' '+cls:'')+'"><div class="ms-l">'+msbEsc(label)+'</div><div class="ms-v">'+msbEsc(vals)+'</div><div class="ms-r">'+msbEsc(res)+'</div><div class="ms-t">'+msbEsc(tag)+'</div></div>';
}
function msbNotVerifiedRows(list){
  return (list||[]).map(m=>'<div class="ms-row ms-nv"><div class="ms-l">NOT VERIFIED</div><div class="ms-v ms-nv-msg">'+msbEsc(m)+'</div><div class="ms-r"></div><div class="ms-t"><span class="ms-warn">NOT VERIFIED</span></div></div>').join('');
}
// 20 Sep 2026 review: advisory rows (engine notes that do not enter the verdict), one per message
function msbAdvisoryRows(list){
  // 20 Sep 2026: an engine note opening with "FAIL:" (k_alpha unbounded, M_y,Ed >= M_cr) is a failure row, printed with the Warning tag
  return (list||[]).map(m=> /^FAIL:/.test(String(m))
    ? '<div class="ms-row ms-nv"><div class="ms-l">FAIL</div><div class="ms-v ms-nv-msg">'+msbEsc(String(m).replace(/^FAIL:\s*/,''))+'</div><div class="ms-r">&ge; 99</div><div class="ms-t"><span class="ms-warn">Warning</span></div></div>'
    : '<div class="ms-row ms-advrow"><div class="ms-l">Advisory</div><div class="ms-v ms-adv-msg">'+msbEsc(m)+'</div><div class="ms-r"></div><div class="ms-t">advisory</div></div>').join('');
}
// 20 Sep 2026 review: the eigen route's buckled mode shape (LT.mode: x mm, phi and v normalised to their own
// peaks; LT.vPoints the lateral restraint stations), formerly a figure of the deleted CED report. Pure SVG string,
// no hover data (the curves are shapes, not values); empty when the route has no mode.
function msbModeShape(LT){
  const md=LT && LT.mode;
  if(!md || !md.x || md.x.length<3) return '';
  const W=540,H=120,pd=12,Lm=md.x[md.x.length-1]||1;
  let vmx=0; md.v.forEach(vv=>{ vmx=Math.max(vmx,Math.abs(vv)); });
  const X=x=>pd+(W-2*pd)*x/Lm, Yc=H/2, ampl=H/2-pd-14;
  let pphi='', pv='';
  md.x.forEach((x,i)=>{ pphi+=(i? ' L ':'M ')+X(x).toFixed(1)+' '+(Yc-ampl*md.phi[i]).toFixed(1); pv+=(i? ' L ':'M ')+X(x).toFixed(1)+' '+(Yc-ampl*(vmx>0? md.v[i]/vmx : 0)).toFixed(1); });
  const marks=(LT.vPoints||[]).map(x=>{ const xx=X(x).toFixed(1); return '<line x1="'+xx+'" y1="'+pd+'" x2="'+xx+'" y2="'+(H-pd)+'" stroke="#b91c1c" stroke-width="1" stroke-dasharray="3,3"/><text x="'+xx+'" y="'+(H-2)+'" font-family="Arial" font-size="9" text-anchor="middle" fill="#b91c1c">'+g(x/1000,2)+'</text>'; }).join('');
  return '<div class="ms-diag-full ms-mode"><div class="ms-dt">Buckled mode shape &mdash; critical eigenmode (normalised; twist &phi; solid, lateral v dashed; red: lateral restraint points; peak twist at x = '+(LT.modePeakX!=null? g(LT.modePeakX/1000,2) : '&mdash;')+' m)</div>'+
    '<svg class="diag" viewBox="0 0 '+W+' '+H+'" xmlns="http://www.w3.org/2000/svg"><line x1="'+pd+'" y1="'+Yc+'" x2="'+(W-pd)+'" y2="'+Yc+'" stroke="#9ca3af" stroke-width="1"/>'+marks+
    '<path d="'+pphi+'" fill="none" stroke="#1d4ed8" stroke-width="2"/><path d="'+pv+'" fill="none" stroke="#059669" stroke-width="1.6" stroke-dasharray="6,4"/></svg></div>';
}
function msbHead(t){ return '<div class="ms-h">'+t+'</div>'; }
function msbSub(t){ return '<div class="ms-sub">'+t+'</div>'; }
// [20 Sep 2026] reference software second title line when torsion is designed: an open
// section names the warping condition of its ends (the P385 closed forms assume
// fork ends free to warp; the warping-torsion FE reads the end W flags), a
// hollow section simply "Includes Design for Torsion".
function msbTorsionTitle(T,sec){
  if(!T) return '';
  if(sec.isBox || T.box) return '<div>Includes Design for Torsion</div>';
  const e1=(S.ends&&S.ends.e1)||{}, e2=(S.ends&&S.ends.e2)||{};
  const w1=!!(T.fe && e1.warp), w2=!!(T.fe && e2.warp);
  const txt=(w1&&w2)? 'End Warping Fixed' : (!w1&&!w2)? 'Ends Free to Warp' : 'End '+(w1?1:2)+' Warping Fixed, End '+(w1?2:1)+' Free to Warp';
  return '<div>Includes Design for Torsion with Span Warping, '+txt+'</div>';
}
/* ---- [20 Sep 2026] diagram panel of the Member Loading block ----
   The loading sketch (beamDiagram) full width, then the shear force, bending
   moment (tension side down), deflection and - with torsion active - torsional
   moment plots in a two-column grid. plot() (js/05-diagrams.js) rides every
   sample on the <svg> as data-* attributes and appends the hidden hover group
   that installDiagramHover() drives, so each diagram reads "x = .. m  M = ..
   kN.m" under the pointer; opt.name labels the readout. Pure: a.diag and
   a.tors.diag only, nothing recomputed. */
function msbDiagramPanel(a){
  if(typeof plot!=='function') return '';
  const cap=(s)=>'<div class="ms-dt">'+s+'</div>';
  const torsOn=!!(a.tors && a.tors.on && a.tors.diag);
  let h='<div class="ms-diagrams">';
  if(typeof beamDiagram==='function') h+='<div class="ms-diag-full">'+cap('Loading')+beamDiagram(a)+'</div>';
  h+='<div class="ms-diag-grid">'+
    '<div>'+cap('Shear force V (kN)')+plot(a.diag.xs,a.diag.V,{color:'#1d4ed8',fill:'#bcd0f7',unit:'kN',name:'V',fmt:v=>f1(v,2)})+'</div>'+
    '<div>'+cap('Bending moment M (kN.m)')+plot(a.diag.xs,a.diag.M,{color:'#b91c1c',fill:'#f3c2c2',unit:'kN.m',name:'M',flip:true,fmt:v=>f1(v,2)})+'</div>'+
    '<div>'+cap('Deflection &delta; (mm, '+msbEsc(a.governD.combo.label)+')')+plot(a.diag.dx,a.diag.dw,{color:'#166534',fill:'#bfe3cb',unit:'mm',name:'\u03b4',fmt:v=>f1(v,2)})+'</div>'+
    (torsOn? '<div>'+cap('Torsional moment T (kN.m, '+msbEsc(a.tors.governT)+')')+plot(a.tors.diag.xs,a.tors.diag.T,{color:'#7a4',fill:'#dcebc4',unit:'kN.m',name:'T',fmt:v=>f1(v,2)})+'</div>' : '')+
    '</div>';
  h+='<div class="ms-note ms-diag-note"><i>Bending-moment diagram drawn on the tension side (sagging down); shear and deflection to true sign (down = below the axis); hover over a diagram for the value at any point.</i></div>';
  return h+'</div>';
}

/* ---- Web Transverse Forces (EN 1993-1-5 cl 6 + 7.2) block, pure ----
   Prints the web geometry and m1/m2, the full derivation of the governing
   station (worst F_Ed/F_Rd), one table row per station, a stiffener-declared
   advisory row per declared station and the assumptions note. Every number
   is a field of c.web (webTransverseCheck); nothing is recomputed. */
function msbWebTypeText(t){ return t==='a'? 'Fig 6.1(a) interior' : t==='b'? 'Fig 6.1(b) through the web' : 'Fig 6.1(c) end'; }
function msbWebBlock(a,c,sec,nvRows){
  const W=c.web||null;
  if(!W) return msbNotVerifiedRows(nvRows);   // 20 Sep 2026: no block when the engine did not run the check (the rows, if any, still print)
  let h=msbHead('Web Transverse Forces (EN 1993-1-5 cl 6)');
  const webs= W.nWebs>1? W.nWebs+' webs' : '1 web';
  h+=msbRow('Web h<sub>w</sub>, t<sub>w</sub>, t<sub>f</sub>, b<sub>f</sub>', msbMM(W.hw)+', '+msbMM(W.tw)+', '+msbMM(W.tf)+', '+msbMM(W.bf)+' mm ('+(W.isBox? 'B/2' : 'B')+' = '+msbMM(W.bfRaw)+' &le; '+(W.isBox||W.chan? 't<sub>w</sub> + 15&epsilon;t<sub>f</sub>' : 't<sub>w</sub> + 30&epsilon;t<sub>f</sub>')+' = '+msbMM(W.bfLim)+'); f<sub>yw</sub> = f<sub>yf</sub> = '+msbInt(W.fyw)+'; '+webs+(W.isBox? ' (flat depth from the section table, corner geometry)' : ''), '', 'Fig 5.1');
  h+=msbRow('m<sub>1</sub> = f<sub>yf</sub>.b<sub>f</sub>/(f<sub>yw</sub>.t<sub>w</sub>) ; m<sub>2</sub> = 0.02(h<sub>w</sub>/t<sub>f</sub>)&sup2;', msbInt(W.fyf)+' x '+msbMM(W.bf)+'/('+msbInt(W.fyw)+' x '+msbMM(W.tw)+') = '+msbR(W.m1)+' ; 0.02 x ('+msbMM(W.hw)+'/'+msbMM(W.tf)+')&sup2; = '+msbR(W.m2full)+' if &lambda;&#772;<sub>F</sub> &gt; 0.5, else 0', '', '6.5(1)');
  h+=msbRow('a = distance between transverse stiffeners', msbEsc(W.aBasis), '', '6.4(1)');
  const G=W.show||W.gov2||null;
  if(G){
    const s=G, t=s.gov, cs=s.cases[s.g2], c72=s.cases[s.g72];
    const nvTag='<span class="ms-warn">NOT VERIFIED</span>';
    h+=msbSub((s.nv? 'Worst station (NOT VERIFIED: s<sub>s</sub> not entered, lower bound 0) x = ' : 'Governing station x = ')+msbM(s.x/1000)+' m: '+msbEsc(s.label)+', load type ('+t.type+') '+msbWebTypeText(t.type)+(s.types.length>1? ' [types '+s.types.map(x=>'('+x+')').join(', ')+' evaluated, lower F<sub>Rd</sub> governs]' : ''));
    const ssTxt='s<sub>s</sub> = '+msbMM(s.ss)+' mm'+(s.ssDefault? ' (not entered: lower bound 0)' : s.kind==='load'&&s.ssIn===0? ' (default 0)' : ' (entered)')+(s.ssCap? ' (capped at h<sub>w</sub>, 6.3(1))' : '')+'; c = '+msbMM(s.c)+' mm (d = '+msbMM(s.d)+' to the member end)';
    const kfTxt= t.type==='c'? 'k<sub>F</sub> = 2 + 6(s<sub>s</sub> + c)/h<sub>w</sub> &le; 6 = 2 + 6 x '+msbMM(s.ss+s.c)+'/'+msbMM(W.hw) : 'k<sub>F</sub> = '+(t.type==='b'? '3.5' : '6')+' + 2(h<sub>w</sub>/a)&sup2; = '+(t.type==='b'? '3.5' : '6')+' + 2('+msbMM(W.hw)+'/'+g(s.a,0)+')&sup2;';
    h+=msbRow('s<sub>s</sub>, c ; k<sub>F</sub>', ssTxt+' ; '+kfTxt, msbR(t.kF), 'Fig 6.1('+t.type+')');
    h+=msbRow('F<sub>cr</sub> = 0.9k<sub>F</sub>.E.t<sub>w</sub>&sup3;/h<sub>w</sub>', '0.9 x '+msbR(t.kF)+' x '+msbInt(a.E)+' x '+msbMM(W.tw)+'&sup3;/'+msbMM(W.hw), msbKN(t.Fcr)+' kN', '6.4(1)');
    const m2Txt='m<sub>2</sub> = '+msbR(t.m2)+(t.iter? ' (first pass &lambda;&#772;<sub>F</sub> = '+msbR(t.lam1)+' &le; 0.5 with m<sub>2</sub> = '+msbR(W.m2full)+', so m<sub>2</sub> = 0)' : ' (&lambda;&#772;<sub>F</sub> &gt; 0.5)');
    if(t.type==='c'){
      h+=msbRow('l<sub>e</sub> = k<sub>F</sub>.E.t<sub>w</sub>&sup2;/(2f<sub>yw</sub>.h<sub>w</sub>) &le; s<sub>s</sub> + c', msbR(t.kF)+' x '+msbInt(a.E)+' x '+msbMM(W.tw)+'&sup2;/(2 x '+msbInt(W.fyw)+' x '+msbMM(W.hw)+') = '+msbMM(t.leRaw)+' &le; '+msbMM(s.ss+s.c), msbMM(t.le)+' mm', '6.5(4)');
      h+=msbRow('l<sub>y</sub> = min[l<sub>e</sub> + t<sub>f</sub>&radic;(m<sub>1</sub>/2 + (l<sub>e</sub>/t<sub>f</sub>)&sup2; + m<sub>2</sub>), l<sub>e</sub> + t<sub>f</sub>&radic;(m<sub>1</sub> + m<sub>2</sub>)]', m2Txt+'; min['+msbMM(t.l1)+', '+msbMM(t.l2)+']', msbMM(t.ly)+' mm', '6.5(4)');
    } else {
      h+=msbRow('l<sub>y</sub> = s<sub>s</sub> + 2t<sub>f</sub>(1 + &radic;(m<sub>1</sub> + m<sub>2</sub>)) &le; a', m2Txt+'; '+msbMM(s.ss)+' + 2 x '+msbMM(W.tf)+' x (1 + &radic;('+msbR(W.m1)+' + '+msbR(t.m2)+')) = '+msbMM(t.l1)+(t.capA? ' &gt; a = '+g(s.a,0)+' (capped)' : ''), msbMM(t.ly)+' mm', '6.5(3)');
    }
    h+=msbRow('&lambda;&#772;<sub>F</sub> = &radic;(l<sub>y</sub>.t<sub>w</sub>.f<sub>yw</sub>/F<sub>cr</sub>)', '&radic;('+msbMM(t.ly)+' x '+msbMM(W.tw)+' x '+msbInt(W.fyw)+'/'+msbKN(t.Fcr)+'e3)', msbR(t.lam), '6.4(1)');
    h+=msbRow('&chi;<sub>F</sub> = 0.5/&lambda;&#772;<sub>F</sub> &le; 1 ; L<sub>eff</sub> = &chi;<sub>F</sub>.l<sub>y</sub>', '0.5/'+msbR(t.lam)+' = '+msbR(t.chiRaw)+(t.chiRaw>1? ' &rarr; 1.000' : '')+' ; '+msbR(t.chi)+' x '+msbMM(t.ly), msbMM(t.Leff)+' mm', '6.4(1)');
    h+=msbRow('F<sub>Rd</sub> = f<sub>yw</sub>.L<sub>eff</sub>.t<sub>w</sub>/&gamma;<sub>M1</sub>', msbInt(W.fyw)+' x '+msbMM(t.Leff)+' x '+msbMM(W.tw)+'/1.0'+(W.isBox? ' per web; load share to this web = '+msbR(s.share)+(s.eMax>0? ' (lever rule, e = '+g(s.eMax,0)+' mm)' : ' (e = 0)')+'; F<sub>Rd</sub> for the load = '+msbKN(s.FRdTot)+' kN' : ''), msbKN(t.FRd)+' kN', '6.2(1)');
    const fTxt=(cs.P!==0&&s.support? 'P = '+msbKN(Math.abs(cs.P))+', R = '+msbKN(cs.R)+': ' : '')+'F<sub>Ed</sub> = '+msbKN(cs.F)+' kN ('+msbEsc(cs.combo)+'; on the '+cs.flange+' flange, '+msbEsc(cs.flangeState)+')';
    h+=msbRow('F<sub>Ed</sub>/F<sub>Rd</sub>', fTxt+' / '+msbKN(s.FRdTot)+' =', msbR(cs.eta2), s.nv? nvTag : msbWarn(cs.eta2<=1.0001));
    const eta1Txt='M<sub>Ed</sub>/M<sub>c.y.Rd</sub> = '+msbKNm(c72.M)+'/'+msbKNm(W.McRd0)+(W.NEd>1e-9? ' + N<sub>Ed</sub>/N<sub>pl.Rd</sub> = '+msbKN(W.NEd)+'/'+msbKN(W.NplRd) : '')+' = '+msbR(c72.eta1);
    h+=msbRow('&eta;<sub>2</sub> + 0.8&eta;<sub>1</sub> &le; 1.4', msbR(c72.eta2)+' + 0.8 x ('+eta1Txt+') = '+msbR(c72.u72raw)+' &le; 1.4 ('+msbEsc(c72.combo)+(c72.flangeComp? '' : '; loaded flange in tension: 7.2(2) refers to 6.2.1(5), expression applied as a screen [verify]')+')', msbR(c72.u72), (s.nv? nvTag : msbWarn(c72.u72<=1.0001))+' 7.2');
  }
  // every station
  if(W.stations.length){
    h+='<table class="ms-combos"><thead><tr><th>x (m)</th><th>Station</th><th>Type</th><th>s<sub>s</sub> (mm)</th><th>k<sub>F</sub></th><th>l<sub>y</sub> (mm)</th><th>&lambda;&#772;<sub>F</sub></th><th>&chi;<sub>F</sub></th><th>F<sub>Rd</sub> (kN)</th><th>F<sub>Ed</sub> (kN)</th><th>Load case</th><th>F<sub>Ed</sub>/F<sub>Rd</sub></th><th>(&eta;<sub>2</sub>+0.8&eta;<sub>1</sub>)/1.4</th><th></th></tr></thead><tbody>'+
      W.stations.map(s=>{
        if(s.stiff) return '<tr><td class="num">'+msbM(s.x/1000)+'</td><td>'+msbEsc(s.label)+'</td><td colspan="11">'+msbEsc(s.msg)+'</td><td>advisory</td></tr>';
        const t=s.gov, cs=s.cases[s.g2];
        return '<tr><td class="num">'+msbM(s.x/1000)+'</td><td>'+msbEsc(s.label)+'</td><td>('+t.type+')</td><td class="num">'+msbMM(s.ss)+(s.ssDefault? '*' : '')+'</td><td class="num">'+msbR(t.kF)+'</td><td class="num">'+msbMM(t.ly)+'</td><td class="num">'+msbR(t.lam)+'</td><td class="num">'+msbR(t.chi)+'</td><td class="num">'+msbKN(s.FRdTot)+'</td><td class="num">'+msbKN(cs.F)+'</td><td>'+msbEsc(cs.combo)+'</td><td class="num">'+msbR(s.eta2)+'</td><td class="num">'+msbR(s.u72)+'</td><td>'+(s.nv? '<span class="ms-warn">NOT VERIFIED</span>' : (s.eta2<=1.0001&&s.u72<=1.0001)? (s===W.gov2? 'governs' : 'OK') : '<span class="ms-warn">Warning</span>')+'</td></tr>';
      }).join('')+'</tbody></table>';
  }
  h+='<div class="ms-note">'+(W.anyDefaultSs? '* s<sub>s</sub> not entered at this end: evaluated at the lower bound s<sub>s</sub> = 0 (F<sub>Rd</sub> rises with the seating length, so a station passing at 0 is verified for any seating; one failing at 0 is NOT VERIFIED until s<sub>s</sub> is entered). ' : '')+'Point loads act on the top flange (bottom flange for an upward load), reactions on the bottom flange; s<sub>s</sub> &le; h<sub>w</sub> (6.3(1)); c = distance from the bearing edge to the member end; type (c) is evaluated whenever s<sub>s</sub> + c &lt; 2h<sub>w</sub>/3 (k<sub>F</sub>(c) &lt; 6) and the lower F<sub>Rd</sub> of types (a) and (c) governs; a load over an end holding U<sub>z</sub> is type (b) with F<sub>Ed</sub> = max(P, R). &eta;<sub>1</sub> uses the unreduced M<sub>c.y.Rd</sub> ('+(W.cls<=2? 'W<sub>pl.y</sub>' : 'W<sub>el.y</sub>')+') [verify: EN 1993-1-5 4.6 writes &eta;<sub>1</sub> with W<sub>eff</sub>]. Distributed loads, hanger loads, the closely-spaced total-load check (6.3(3)) and flange-induced buckling (section 8) are not evaluated.</div>';
  h+=msbNotVerifiedRows(nvRows);
  return h;
}

function renderDesignBrief(a,c,sec){
  sec=sec||a.sec;
  const LT=c.ltb||null, AX=c.ax||null, B=c.buck||null, T=c.tor||null, AN=c.annex||null;
  const fullRest=(S.restraint||'full')==='full';
  const isCant=isCantilever(S);
  const eigen=!!(LT && LT.eigen);
  const ltbChecked=!fullRest && !(LT && LT.na);
  const utils=c.utils||[];
  const failed=utils.some(u=>!Number.isFinite(u.val)||u.val>1.0001);
  const unsupported=c.unsupported||[];
  const verdict=c.pass? 'PASS' : failed? 'FAIL' : 'NOT VERIFIED';
  const titleSuffix=c.pass? '' : ' ('+verdict+')';
  const famLabel = sec.isBox? (S.family==='rhs'? 'RHS ['+(sec.boxType==='CF'?'Cold-formed':'Hot-finished')+']' : 'SHS ['+(sec.boxType==='CF'?'Cold-formed':'Hot-finished')+']') : S.family==='ub'? 'UB' : S.family==='uc'? 'UC' : 'PFC';
  const secStr=msbSecName(sec.key)+' '+famLabel+' ['+S.grade+']';
  const memberName=(S.memberName&&String(S.memberName).trim())? String(S.memberName).trim() : secStr;
  const cls=c.cl.cls;
  const Wy_cm3=c.Wy/1e3;
  const gfb=a.governM.fb;
  const [xa,xb]=msbPortion(a,LT);
  const wholePortion = xa<=1e-6 && Math.abs(xb-a.L)<=1e-6;
  const nULS=a.ulsResults.length;
  const caseM=msbCaseLabel(a.governM.combo,false,a), caseD=msbCaseLabel(a.governD.combo,true,a);
  // blocking messages sorted into their blocks
  const nv={forces:[],class:[],local:[],web:[],compression:[],ltb:[],torsion:[],general:[]};
  unsupported.forEach(m=>{ let k=msbBlockFor(m); if(k==='compression' && !(B && B.Fc>1e-9)) k='local'; if(k==='torsion' && !T) k= /Annex A|k_alpha/.test(m)? 'ltb' : 'general'; nv[k].push(m); });

  let h='';
  /* ---- 5.0 Title ---- */
  h+='<div class="ms-title"><div>'+(AX? 'Axial with Moments (Member)' : 'Beam &amp; Beam-Portion (Member)')+titleSuffix+'</div>'+
     msbTorsionTitle(T,sec)+   // 20 Sep 2026: "Includes Design for Torsion ..." when torsion is active
     '<div>Member '+memberName+'</div>'+
     '<div>'+(isCant? 'Cantilever 0 to '+msbM(a.L/1000)+' m' : 'Between '+msbM(xa/1000)+' and '+msbM(xb/1000)+' m')+', in Load Case '+caseM+'</div>'+
     // [beam-v03 addition, 19 Sep 2026 scope] the seven DOF flags of each end, the derived end types, the preset name and the hinges
     '<div class="ms-ends">End conditions: '+endsConditionsLine(S)+'</div></div>';

  /* ---- 5.1 Member Loading and Member Forces ---- */
  h+=msbHead('Member Loading and Member Forces');
  const loadLines=[];
  loadLines.push('Loading Combination : '+a.governM.combo.label+(nULS>1? ' (governing moment; '+nULS+' ULS cases enabled)' : ''));
  const ecc=(ld)=> S.eccOn && ld.type!=='moment'? ' e = '+g(ld.e||0,0)+' mm'+(loadHeightPerLoadOn()? ' z<sub>g</sub> = '+g(loadZgValue(ld),0)+' mm' : '') : '';
  const lab=(ld)=> ld.label? ' &mdash; '+String(ld.label).replace(/[&<>]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[m])) : '';   // 20 Sep 2026: optional load label / origin text (wall.html ledgers write "Inner leaf: floor joists")
  S.loads.filter(ld=>!ld.isSelfWeight).forEach(ld=>{
    const cs=ld.case||'?';
    if(ld.type==='udl') loadLines.push(cs+' UDL '+f1(ld.w,3)+' '+g(ld.x1)+'&ndash;'+g(ld.x2)+' m'+ecc(ld)+' ( kN/m )'+lab(ld));
    else if(ld.type==='trap') loadLines.push(cs+' TRAP '+f1(ld.w1,3)+'&rarr;'+f1(ld.w2,3)+' '+g(ld.x1)+'&ndash;'+g(ld.x2)+' m'+ecc(ld)+' ( kN/m )'+lab(ld));
    else if(ld.type==='point') loadLines.push(cs+' PY '+f1(ld.P,3)+' @ '+g(ld.pos)+' m'+ecc(ld)+' ( kN )'+lab(ld));
    else if(ld.type==='moment') loadLines.push(cs+' M '+f1(ld.M,3)+' @ '+g(ld.pos)+' m ( kN.m )'+lab(ld));
  });
  const swE=selfWeightEccentricity(sec);
  loadLines.push('G SW '+g(selfWeightValue(sec),4)+' kN/m 0&ndash;'+g(S.L)+' m'+(S.eccOn&&Math.abs(swE)>1e-9? ' e = '+g(swE,1)+' mm' : '')+' ( automatic )');
  // gamma_G,inf companions (19 Sep 2026 review): G at 1.0 (STR set B) and 0.9 (EQU set A) of every ULS combination
  // with G > 1.0, solved for the end reactions (uplift / hold-down, web bearing)
  if(a.ulsCompanions && a.ulsCompanions.length){
    loadLines.push('<b>&gamma;<sub>G,inf</sub> companions</b> (reactions only): '+a.ulsCompanions.length+' &mdash; G at 1.0 (STR set B) and 0.9 (EQU set A) of every ULS combination with G &gt; 1.0, variable factors as entered');
    a.ulsCompanions.forEach((r,i)=>loadLines.push('&nbsp;&nbsp;ULS C'+(i+1)+': '+r.combo.label));
  }
  if(a.companionNote) loadLines.push('<span class="ms-note">'+a.companionNote+'</span>');
  // 20 Sep 2026: the load list stands alone (the section load-lines figure of
  // js/05-section-view.js beside it when eccentricities / load heights are
  // entered); the sketch and the value diagrams follow the forces table.
  const secFig=(typeof sectionLoadLineView==='function' && (S.eccOn || (typeof loadHeightPerLoadOn==='function' && loadHeightPerLoadOn())))? sectionLoadLineView(sec) : '';
  h+='<div class="ms-loading'+(secFig? '' : ' ms-loading-plain')+'"><div class="ms-loadlist">'+loadLines.join('<br>')+'</div>'+(secFig? '<div class="ms-sketch">'+secFig+'</div>' : '')+'</div>';

  // member forces table
  const N=S.axial||0, Ntag=N>=0? 'C':'T';
  // 20 Sep 2026 review: the Torque Moment column is the end torque REACTION, as reference software prints it: the total
  // T = GI_T phi' - EI_w phi''' (T.TEnds) when the P385 / FE torsion analysis provides it, else the St Venant part
  const torqueAt=(x)=>{ if(T && T.p385 && T.TEnds) return x<=xa+1e-6? T.TEnds[0] : T.TEnds[1]; if(T && T.p385 && T.TtEnds) return x<=xa+1e-6? T.TtEnds[0] : T.TtEnds[1]; if(a.tors && a.tors.diag) return interpAt(a.tors.diag.xs,a.tors.diag.T,x/1000); return 0; };
  const V1=msbEndShear(gfb,xa+1e-4), V2=msbEndShear(gfb,xb-1e-4);
  const M1=wholePortion? a.M0end : msbMomentAt(gfb,xa+1e-4), M2=wholePortion? a.MLend : msbMomentAt(gfb,xb-1e-4);
  const pm=msbPortionMoments(gfb,xa,xb);
  const MmaxP=wholePortion? a.Mmax : pm.Mmax, MposP=wholePortion? a.Mpos : pm.xmax/1000;
  const dfl=a.deflection||{dmax:a.dmax,dpos:a.dpos*1000};
  // [beam-v03 addition, 19 Sep 2026 scope] End 1 / End 2 reactions of the governing-moment
  // combination (a.reactions: end, type, V, M): R and M as the end type carries them, a
  // guided end "M only", a free end none; printed only when the portion end is the member end
  const reactCells=(n,xEnd)=>{
    const r=a.reactions.find(q=>q.end===n);
    if(Math.abs(xEnd-(n===1? 0 : a.L))>1e-6) return '<td class="num">&mdash;</td><td class="num">&mdash;</td>';
    if(!r) return '<td class="num">free: none</td><td class="num">&mdash;</td>';
    if(r.type==='guided') return '<td class="num">guided: M only</td><td class="num">'+f1(reactionEndMomentKNm(r),2)+'</td>';
    return '<td class="num">'+f1(r.V/1000,2)+'</td><td class="num">'+(r.type==='fixed'? f1(reactionEndMomentKNm(r),2) : '&mdash;')+'</td>';
  };
  h+='<table class="ms-forces"><thead><tr><th colspan="13" class="ms-ft">Member Forces in Load Case '+caseM+' and Maximum Deflection from Load Case '+caseD+'</th></tr>'+
     '<tr><th rowspan="2">Mem<br>ber<br>No.</th><th rowspan="2">Node<br>End1<br>End2</th><th rowspan="2">Axial<br>Force<br>(kN)</th><th rowspan="2">Torque<br>Moment<br>(kN.m)</th><th colspan="2">Shear Force<br>(kN)</th><th colspan="2">Bending Moment<br>(kN.m)</th><th colspan="2">Maximum Moment<br>(kN.m @ m)</th><th rowspan="2">Maximum<br>Deflection<br>(mm @ m)</th><th colspan="2">Reaction<br>R (kN), M (kN.m)</th></tr>'+
     '<tr><th>y-y</th><th>z-z</th><th>y-y</th><th>z-z</th><th>y-y</th><th>z-z</th><th>R</th><th>M</th></tr></thead><tbody>'+
     '<tr><td class="num">1</td><td class="num">End 1, x = '+msbM(xa/1000)+'</td><td class="num">'+f1(Math.abs(N),2)+Ntag+'</td><td class="num">'+f1(torqueAt(xa),2)+'</td><td class="num">'+f1(V1,2)+'</td><td class="num">0.00</td><td class="num">'+f1(M1,2)+'</td><td class="num">'+f1(S.Mz||0,2)+'</td><td class="num">'+f1(MmaxP,2)+'</td><td class="num">'+f1(S.Mz||0,2)+'</td><td class="num">'+f1(dfl.dmax,2)+'</td>'+reactCells(1,xa)+'</tr>'+
     '<tr><td class="num"></td><td class="num">End 2, x = '+msbM(xb/1000)+'</td><td class="num">'+f1(Math.abs(N),2)+Ntag+'</td><td class="num">'+f1(torqueAt(xb),2)+'</td><td class="num">'+f1(V2,2)+'</td><td class="num">0.00</td><td class="num">'+f1(M2,2)+'</td><td class="num">'+f1(S.Mz||0,2)+'</td><td class="num">@ '+f1(MposP,3)+'</td><td class="num">@ &mdash;</td><td class="num">@ '+f1(dfl.dpos/1000,3)+'</td>'+reactCells(2,xb)+'</tr>'+
     '</tbody></table>';
  h+='<div class="ms-note">Reactions of '+a.governM.combo.label+' at the member ends: R positive upward (a negative R is uplift), M = the end bending moment in the diagram convention (sagging positive, hogging negative); a pinned end carries no M, a guided end no R, a free end neither. V<sub>z</sub> = 0: no minor-axis shear in the single-plane model; M<sub>z</sub> is the entered constant design moment.</div>';
  // uplift / hold-down (item 1.2): one row per lifting support, in the Member Forces block
  const HD=c.holdDown||null;
  if(HD && HD.rows && HD.rows.length){
    HD.rows.forEach(u=>{
      if(u.level==='sls'){
        h+=msbRow('Uplift at SLS only, End '+u.n+' (x = '+msbM(u.pos/1000)+' m)', 'R = &minus;'+f1(Math.abs(u.RSls),2)+' kN (SLS combination '+msbEsc(u.comboSls)+'); no ULS combination lifts this end, incl. the &gamma;<sub>G,inf</sub> companions (G at 1.0 STR set B, 0.9 EQU set A) (advisory)', 'R = &minus;'+f1(Math.abs(u.RSls),2)+' kN', 'advisory');
        return;
      }
      const lbl='Hold-down '+(u.holdDown? 'provided' : 'required')+' at End '+u.n+' (x = '+msbM(u.pos/1000)+' m)';
      const vals='R = &minus;'+f1(Math.abs(u.RUls),2)+' kN (combination '+msbEsc(u.comboUls)+')'+(u.RSls!=null? '; SLS uplift &minus;'+f1(Math.abs(u.RSls),2)+' kN ('+msbEsc(u.comboSls)+')' : '')+'; '+u.nCombos+' combination(s) lift this end';
      if(u.holdDown) h+=msbRow(lbl, vals+' &mdash; design the hold-down connection for this force (advisory)', 'R = &minus;'+f1(Math.abs(u.RUls),2)+' kN', 'hold-down');
    });
    h+=msbNotVerifiedRows(nv.forces);
  }
  if(a.uplift && !a.uplift.any) h+=msbRow('Uplift', 'no end lifts in any of the '+(a.uplift.nCombos!=null? a.uplift.nCombos : a.ulsResults.length+a.slsResults.length)+' combinations (all reactions &ge; 0'+((a.ulsCompanions&&a.ulsCompanions.length)? '; incl. the '+a.ulsCompanions.length+' &gamma;<sub>G,inf</sub> companions with G at 1.0 and 0.9' : '')+')', 'R<sub>min</sub> &ge; 0', 'OK');
  // 20 Sep 2026: loading sketch + hover diagrams (V, M, delta, T) after the forces table
  h+=msbDiagramPanel(a);
  if(nULS>1 || a.slsResults.length>1){
    h+='<table class="ms-combos"><thead><tr><th>Combination</th><th>V<sub>max</sub> (kN)</th><th>M<sub>max</sub> (kN.m @ m)</th><th>&delta; (mm)</th></tr></thead><tbody>'+
       a.ulsResults.map(r=>'<tr><td>'+r.combo.label+(r===a.governM? ' (governs M)':'')+'</td><td class="num">'+f1(r.Vmax/1000,3)+'</td><td class="num">'+f1(r.Mmax/1e6,3)+' @ '+g(r.Mpos/1000,3)+'</td><td class="num">&mdash;</td></tr>').join('')+
       a.slsResults.map(r=>'<tr><td>'+r.combo.label+(r===a.governD? ' (governs &delta;)':'')+'</td><td class="num">&mdash;</td><td class="num">&mdash;</td><td class="num">'+f1(r.dmax,2)+'</td></tr>').join('')+
       '</tbody></table>';
  }

  /* ---- 5.2 Classification and Effective Area ---- */
  h+=msbHead('Classification and Effective Area (EN 1993: 2006)');
  h+=msbRow('Section ('+f1(sec.mass,1)+' kg/m)', secStr+(sec.isBox? ' D='+g(sec.D,0)+' B='+g(sec.B,0)+' t='+g(sec.tf,1) : ''),'','');
  const clsTag= cls>=4? '<span class="ms-warn">BLOCKED</span>' : 'Class '+cls;
  h+=msbRow('Class = Fn(c/t, d/t, f<sub>y</sub>, N, M<sub>y</sub>, M<sub>z</sub>)',
    f1(sec.bT,2)+', '+f1(sec.dt,2)+', '+msbInt(c.fy)+', '+f1(Math.abs(N),2)+', '+f1(c.Mx,2)+', '+f1(Math.abs(S.Mz||0),2)+'<span class="ms-inline-tag">'+msbAxialTag(sec,c.eps)+'</span>',
    '', clsTag);
  if((c.cl.webCase && c.cl.webCase!=='bending') || c.cl.mzStress){
    const wl=c.cl.wlim||[];
    const mzWeb = c.cl.mzStress? ' (M<sub>z</sub> does not stress the web of an I/H section: it lies on the z-z axis)' : '';
    h+=msbRow('Web classified for', c.cl.webCase==='bending+compression'
      ? 'bending + compression: &alpha; = '+f1(c.cl.alphaW,3)+', &psi; = '+f1(c.cl.psiW,3)+'; limits '+g(wl[0],1)+'&epsilon; / '+g(wl[1],1)+'&epsilon; / '+g(wl[2],1)+'&epsilon;'+mzWeb
      : c.cl.webCase==='bending'? 'y-y bending: limits 72&epsilon; / 83&epsilon; / 124&epsilon;'+mzWeb
      : 'biaxial: uniform-compression web bound; limits '+g(wl[0],1)+'&epsilon; / '+g(wl[1],1)+'&epsilon; / '+g(wl[2],1)+'&epsilon;',
      'Class '+c.cl.wc, 'Table 5.2');
  }
  // [beam-v03 addition, 19 Sep 2026 G3] I/H flange outstands under M_y + M_z: the outstand compressed by M_z is wholly in compression (alpha = 1)
  if(c.cl.mzStress){
    const ms=c.cl.mzStress;
    h+=msbRow('Flange outstands under M<sub>y</sub> + M<sub>z</sub> (+N)', '&sigma; = N/A + M<sub>y</sub>/W<sub>el.y</sub> &plusmn; M<sub>z</sub>y/I<sub>z</sub>: outstand compressed by M<sub>z</sub>: root '+f1(ms.compRoot,1)+', tip '+f1(ms.compTip,1)+' N/mm&sup2; (both compressive, &alpha; = 1 &rarr; 9&epsilon;/&alpha;, 10&epsilon;/&alpha;, 21&epsilon;&radic;k<sub>&sigma;</sub> &ge; 9&epsilon;/10&epsilon;/14&epsilon; bound kept); opposite outstand: root '+f1(ms.relRoot,1)+', tip '+f1(ms.relTip,1)+' N/mm&sup2; ('+ms.relState+', laxer limits, does not govern)', 'Class '+c.cl.fc, 'Table 5.2 sheet 2');
  }
  // [beam-v03 addition, 19 Sep 2026 G3] effective area of a Class-4 web in uniform compression (EN 1993-1-5 4.4)
  if(c.aeff && c.aeff.active){
    const ae=c.aeff;
    h+=msbRow('A<sub>eff</sub> = A &minus; '+(ae.nWebs>1? '2':'')+'(1 &minus; &rho;)b&#772;t<sub>w</sub>', 'web Class 4 in uniform compression: d/t = '+f1(ae.dt,2)+' &gt; 42&epsilon; = '+f1(ae.limit,2)+'; &lambda;&#772;<sub>p</sub> = (b&#772;/t)/(28.4&epsilon;&radic;k<sub>&sigma;</sub>) = '+f1(ae.dt,2)+'/(28.4 x '+f1(ae.eps,3)+' x &radic;4) = '+msbR(ae.lamP)+'; &rho; = (&lambda;&#772;<sub>p</sub> &minus; 0.055(3+&psi;))/&lambda;&#772;<sub>p</sub>&sup2; = '+msbR(ae.rho)+' (&psi; = 1); b<sub>eff</sub> = '+msbMM(ae.beff)+' mm = b<sub>e1</sub> + b<sub>e2</sub> = 2 x '+msbMM(ae.be1)+'; '+f1(ae.A/100,2)+' &minus; '+(ae.nWebs>1? '2 x ':'')+f1(ae.bineff*ae.tw/100,2)+' cm&sup2; (e<sub>N</sub> = 0, symmetric)', f1(ae.Aeff/100,2)+' cm&sup2;', 'EN 1993-1-5 4.4');
  }
  const ulsIdx=a.ulsResults.map((r,i)=>i+1), slsIdx=a.slsResults.map((r,i)=>i+1);
  h+=msbRow('Auto Design Load Cases', msbCaseRanges(ulsIdx)+(slsIdx.length? '; SLS '+msbCaseRanges(slsIdx) : ''),'','');
  h+=msbNotVerifiedRows(nv.class);

  /* ---- 5.2b Shear Capacity Check (reference software block; 20 Sep 2026 order) ----
     the maximum shear of every ULS combination against V_pl.y.Rd (cl 6.2.6)
     and the cl 6.2.6(6) shear-buckling screen; no V_z.Ed line because the
     single-plane model has no minor-axis shear (reference software prints one only
     with a minor-axis load) */
  h+=msbHead('Shear Capacity Check');
  h+=msbRow('V<sub>pl.y.Rd</sub> = A<sub>v</sub>f<sub>y</sub>/(&radic;3&gamma;<sub>M0</sub>)', g(c.Av,1)+' mm&sup2; x '+msbInt(c.fy)+'/(&radic;3 x 1)'+(c.avFloor!=null? ' ; A<sub>v</sub> &ge; &eta;h<sub>w</sub>t<sub>w</sub> = '+g(c.avFloor,1)+' mm&sup2;' : ''), msbKN(c.VcRd)+' kN', '6.2.6');
  h+=msbRow('V<sub>y.Ed</sub>/V<sub>pl.y.Rd</sub>', msbKN(c.Fv)+' / '+msbKN(c.VcRd)+' =', msbR(c.shearUtil), msbWarn(c.shearUtil<=1.0001));
  if(c.sbRatio!=null && c.sbLimit!=null) h+=msbRow((sec.isBox? 'd/t' : 'h<sub>w</sub>/t<sub>w</sub>')+' &le; 72&epsilon;/&eta;', f1(c.sbRatio,2)+' '+(c.sbOk? '&le;' : '&gt;')+' 72 x '+f1(c.eps,3)+'/'+g(c.eta!=null? c.eta : 1,2)+' = '+f1(c.sbLimit,2)+(c.sbOk? '' : ' (shear buckling, EN 1993-1-5 5: not evaluated)'), c.sbOk? 'no shear buckling' : '&mdash;', c.sbOk? 'OK 6.2.6(6)' : '<span class="ms-warn">BLOCKED</span>');

  /* ---- 5.3 Local Capacity Check / Moment Capacity Check ---- */
  h+=msbHead(AX? 'Local Capacity Check' : 'Moment Capacity Check M.c.y.Rd'+(fullRest? ' - Fully Restrained Beam' : ''));
  const VplMoment=(T && T.VplTRd!=null)? T.VplTRd : c.VcRd;
  // the shear coincident with the maximum moment (cl 6.2.8 low-shear test), as reference software prints it here
  h+=msbRow('V<sub>y.Ed</sub>/V<sub>pl.y.Rd</sub> (at max M)', msbKN(c.VatM)+' / '+msbKN(c.VcRd)+' =', msbR(c.VatM/Math.max(c.VcRd,1e-9)), c.lowShearAtM? 'Low Shear' : 'High Shear');   // 20 Sep 2026 review: "(at max M)" in the label, as reference software
  const hsReduced = !c.lowShearAtM && /6\.2\.8\(3\)/.test(msbEsc(c.hsNote));
  if(hsReduced){
    h+=msbRow('&rho; = (2V<sub>y.Ed</sub>/V<sub>pl'+(T&&T.VplTRd!=null?'.T':'')+'.Rd</sub> &minus; 1)&sup2;', '(2 x '+msbKN(c.VatM)+'/'+msbKN(VplMoment)+' &minus; 1)&sup2;', msbR(msbRhoShear(c.VatM,VplMoment)), '6.2.8(3)');
  }
  const Wlbl= cls<=2? 'W<sub>pl.y</sub>' : 'W<sub>el.y</sub>';
  h+=msbRow(hsReduced? 'M<sub>v.y.Rd</sub> = '+(c.mvForm||'(W<sub>pl.y</sub> &minus; &rho;A<sub>v</sub>&sup2;/4t<sub>w</sub>)f<sub>y</sub>/&gamma;<sub>M0</sub>') : 'M<sub>c.y.Rd</sub> = f<sub>y</sub>.'+Wlbl+'/&gamma;<sub>M0</sub>',
    msbInt(c.fy)+' x '+f1(Wy_cm3,1)+'/1'+(hsReduced&&c.rhoAtM!=null? ' with &rho; = '+msbR(c.rhoAtM) : ''), msbKNm(c.McRd)+' kN.m', hsReduced? '6.2.8' : '');
  h+=msbRow('M<sub>y.Ed</sub>/M<sub>c.y.Rd</sub>', msbKNm(c.Mx)+' / '+msbKNm(c.McRd)+' =', msbR(c.momUtil), msbWarn(c.momUtil<=1.0001));
  if(c.coex){
    const cx=c.coex;
    h+=msbRow('M<sub>y.Ed</sub>/M<sub>v.y.Rd</sub> @ x', cx.pureShearFail
      ? '@ x = '+msbM(cx.x/1000)+' m: V<sub>Ed</sub> = '+msbKN(cx.V)+' &gt; V<sub>pl.Rd</sub> = '+msbKN(cx.VplRd)+': pure shear governs'
      : '@ x = '+msbM(cx.x/1000)+' m: M = '+msbKNm(cx.M)+', V = '+msbKN(cx.V)+' &gt; 0.5V<sub>pl.Rd</sub>; &rho; = '+msbR(cx.rho||0)+'; M<sub>v.y.Rd</sub> = '+(cx.form||'')+' = '+msbKNm(cx.MvRd),
      msbR(cx.u), msbWarn(cx.u<=1.0001));
  }
  // [beam-v03 addition, 19 Sep 2026 G3] cl 6.2.10: bending + shear + axial / minor-axis at the worst high-shear station
  if(c.mvn){
    const m=c.mvn;
    const lbl= m.plastic? (m.biax? '(M<sub>y.Ed</sub>/M<sub>N.V.y.Rd</sub>)<sup>&alpha;</sup>+(M<sub>z.Ed</sub>/M<sub>N.V.z.Rd</sub>)<sup>&beta;</sup> @ x' : 'M<sub>y.Ed</sub>/M<sub>N.V.y.Rd</sub> @ x')
                        : 'N<sub>Ed</sub>/N<sub>V.Rd</sub> + M<sub>y.Ed</sub>/M<sub>v.y.Rd</sub>'+(m.biax? ' + M<sub>z.Ed</sub>/M<sub>v.z.Rd</sub>' : '')+' @ x';
    const vals='@ x = '+msbM(m.x/1000)+' m ('+msbEsc(m.combo)+'): V = '+msbKN(m.V)+', &rho; = '+msbR(m.rho)+'; N<sub>V.Rd</sub> = (A &minus; &rho;A<sub>v</sub>)f<sub>y</sub> = '+msbKN(m.NV)+', n<sub>V</sub> = '+msbR(m.nV)+'; M<sub>v.y.Rd</sub> = '+msbEsc(m.formY)+' = '+msbKNm(m.MvY)+(m.biax? '; M<sub>v.z.Rd</sub> = '+msbEsc(m.formZ)+' = '+msbKNm(m.MvZ) : '')+
      (m.plastic? '; a<sub>V</sub> = '+msbR(m.aV)+(m.afV!=null? ', a<sub>f.V</sub> = '+msbR(m.afV) : '')+(m.waiver? '; 6.2.9.1(4): no reduction (N &le; 0.25N<sub>V.Rd</sub>, N &le; 0.5h<sub>w</sub>t<sub>w</sub>(1&minus;&rho;)f<sub>y</sub>)' : '')+'; M<sub>N.V.y.Rd</sub> = '+msbKNm(m.MNVy)+(m.biax? ', M<sub>N.V.z.Rd</sub> = '+msbKNm(m.MNVz)+', &alpha; = '+g(m.alpha,2)+', &beta; = '+g(m.beta,2) : '') : '')+
      '; M<sub>y</sub> = '+msbKNm(m.M)+(m.biax? ', M<sub>z</sub> = '+msbKNm(m.Mz) : '')+', N = '+msbKN(m.N)+' ('+msbEsc(m.form)+')';
    h+=msbRow(lbl, vals, msbR(m.u), msbWarn(m.u<=1.0001)+' 6.2.10');
  }
  if(AX){
    // V_pl.z.Rd, A_v,z and M_c.z.Rd are engine values (c.ax, cl 6.2.6(3) / 6.2.5);
    // 20 Sep 2026: printed with a minor-axis moment only, as reference software does
    if(AX.biax){
      if(AX.VplZ!=null) h+=msbRow('V<sub>z.Ed</sub>/V<sub>pl.z.Rd</sub>', '0 / '+msbKN(AX.VplZ)+' = (A<sub>v,z</sub> = '+g(AX.Avz,1)+' mm&sup2;; no minor-axis shear in the single-plane model)', '0.000', 'Low Shear');
      else h+=msbRow('V<sub>z.Ed</sub>/V<sub>pl.z.Rd</sub>', 'V<sub>pl.z.Rd</sub> not evaluated by the engine for this section', '&mdash;', 'not evaluated');
      if(AX.Mcz!=null) h+=msbRow('M<sub>c.z.Rd</sub> = f<sub>y</sub>.'+(cls<=2? 'W<sub>pl.z</sub>':'W<sub>el.z</sub>')+'/&gamma;<sub>M0</sub>', msbInt(c.fy)+' x '+f1(cls<=2? sec.Sy : sec.Zy,1)+'/1', msbKNm(AX.Mcz)+' kN.m', '');
      else h+=msbRow('M<sub>c.z.Rd</sub>', 'not evaluated by the engine', '&mdash;', 'not evaluated');
    }
    h+=msbRow('N<sub>pl.Rd</sub> = A<sub>g</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(sec.A,2)+' x '+msbInt(c.fy)+'/1 = (No bearing / block tearing design)', msbKN(AX.NplRd)+' kN', '');
    if(AX.aeff && AX.aeff.active) h+=msbRow('N<sub>c.Rd</sub> = A<sub>eff</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(AX.aeff.Aeff/100,2)+' x '+msbInt(c.fy)+'/1 (Class-4 web in uniform compression, cl 6.2.4(2))', msbKN(AX.NcRd)+' kN', '6.2.4');
    if(AX.tension && S.anet!=null) h+=msbRow('N<sub>u.Rd</sub> = 0.9A<sub>net</sub>f<sub>u</sub>/&gamma;<sub>M2</sub>', '0.9 x '+f1(S.anet,2)+' x '+msbInt(fuFromGrade(S.grade))+'/1.10', msbKN(AX.NuRd)+' kN', '&gamma;<sub>M2</sub> = 1.10 (UK NA)');
    if(AX.tension) h+=msbRow('n = N<sub>Ed</sub>/N<sub>t.Rd</sub>', f1(-Math.abs(N),3)+' / '+msbKN(AX.NtRd)+' =', msbR(AX.nUtil), msbWarn(AX.nUtil<=1.0001));
    else if(AX.aeff && AX.aeff.active){
      h+=msbRow('n = N<sub>Ed</sub>/N<sub>pl.Rd</sub>', f1(N,3)+' / '+msbKN(AX.NplRd)+' = (parameter of the 6.2.9 interaction)', msbR(AX.n), '');
      h+=msbRow('N<sub>Ed</sub>/N<sub>c.Rd</sub>', f1(N,3)+' / '+msbKN(AX.NcRd)+' =', msbR(AX.nUtil), msbWarn(AX.nUtil<=1.0001));
    } else h+=msbRow('n = N<sub>Ed</sub>/N<sub>pl.Rd</sub>', f1(N,3)+' / '+msbKN(AX.NplRd)+' =', msbR(AX.n), msbWarn(AX.nUtil<=1.0001));
    if(AX.cls3){
      h+=msbRow('W<sub>el.y</sub>', 'Class 3: elastic (cl 6.2.9.2)', f1(sec.Zx,1)+' cm&sup3;', '');
      h+=msbRow('M<sub>el.y.Rd</sub> = W<sub>el.y</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(sec.Zx,1)+' x '+msbInt(c.fy)+'/1', msbKNm(AX.MN)+' kN.m', '');
      if(AX.biax) h+=msbRow('M<sub>el.z.Rd</sub> = W<sub>el.z</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(sec.Zy,1)+' x '+msbInt(c.fy)+'/1', msbKNm(AX.MNz)+' kN.m', '');
      h+=msbRow('N<sub>Ed</sub>/N<sub>pl.Rd</sub> + M<sub>y.Ed</sub>/M<sub>el.y.Rd</sub> + M<sub>z.Ed</sub>/M<sub>el.z.Rd</sub>', msbR(AX.n)+' + '+msbKNm(c.Mx)+'/'+msbKNm(AX.MN)+' + '+(AX.biax? msbKNm(AX.Mz)+'/'+msbKNm(AX.MNz) : '0')+' =', msbR(AX.mUtil), msbWarn(AX.mUtil<=1.0001));
    } else if(AX.chan){
      h+=msbRow('M<sub>c.y.Rd</sub>, M<sub>c.z.Rd</sub> (channel, no M<sub>N</sub> reduction)', 'cl 6.2.1(7) linear interaction: '+AX.mnLbl, msbKNm(AX.MN)+', '+msbKNm(AX.MNz)+' kN.m', '6.2.1(7)');
      h+=msbRow('n + M<sub>y.Ed</sub>/M<sub>c.y.Rd</sub> + M<sub>z.Ed</sub>/M<sub>c.z.Rd</sub>', msbR(AX.n)+' + '+msbKNm(c.Mx)+'/'+msbKNm(AX.MN)+' + '+(AX.biax? msbKNm(AX.Mz)+'/'+msbKNm(AX.MNz) : '0')+' =', msbR(AX.mUtil), msbWarn(AX.mUtil<=1.0001));
    } else {
      const waiver=/small axial/.test(msbEsc(AX.mnLbl));
      h+=msbRow('W<sub>pl.N.y</sub> = Fn(W<sub>pl.y</sub>, A<sub>vy</sub>, n)', f1(sec.Sx,1)+', '+f1(c.Av/100,3)+', '+msbR(AX.n), f1(msbWplN(AX.MN,c.fy),1)+' cm&sup3;', '');
      h+=msbRow('M<sub>N.y.Rd</sub> = W<sub>pl.N.y</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(msbWplN(AX.MN,c.fy),1)+' x '+msbInt(c.fy)+'/1', msbKNm(AX.MN)+' kN.m', waiver? '6.2.9.1(4)' : '');
      if(AX.biax){
        h+=msbRow('W<sub>pl.N.z</sub> = Fn(W<sub>pl.z</sub>, A<sub>vz</sub>, n)', f1(sec.Sy,1)+', '+(AX.Avz!=null? f1(AX.Avz/100,3) : '&mdash;')+', '+msbR(AX.n), f1(msbWplN(AX.MNz,c.fy),1)+' cm&sup3;', '');
        h+=msbRow('M<sub>N.z.Rd</sub> = W<sub>pl.N.z</sub>.f<sub>y</sub>/&gamma;<sub>M0</sub>', f1(msbWplN(AX.MNz,c.fy),1)+' x '+msbInt(c.fy)+'/1', msbKNm(AX.MNz)+' kN.m', '');
      }
      h+=msbRow('(M<sub>y.Ed</sub>/M<sub>N.y.Rd</sub>)<sup>&alpha;</sup>+(M<sub>z.Ed</sub>/M<sub>N.z.Rd</sub>)<sup>&beta;</sup>',
        '('+msbKNm(c.Mx)+'/'+msbKNm(AX.MN)+')<sup>'+g(AX.alpha,2)+'</sup>+'+(AX.biax? '('+msbKNm(AX.Mz)+'/'+msbKNm(AX.MNz)+')<sup>'+g(AX.beta,2)+'</sup>' : '(0)<sup>1</sup>')+'=',
        msbR(AX.mUtil), msbWarn(AX.mUtil<=1.0001));
    }
  }
  // [20 Sep 2026 torsion + N/Mz] EN 1993-1-1 6.2.7(5): the elastic yield criterion 6.2.1(5) Eq (6.1) with every stress
  // of the torsion solution (c.tor.elastic) - a code check in every torsion case (N_Ed = 0, M_z = 0 included), printed
  // after the plastic N-M-Mz interaction rows; the twist-induced phi.M_y enters M_z,tot (msbElasticRows above)
  if(T && T.elastic) h+=msbElasticRows(T);
  h+=msbNotVerifiedRows(nv.local);

  /* ---- [beam-v03 addition, 19 Sep 2026] Web Transverse Forces (EN 1993-1-5 cl 6), after Local Capacity;
     20 Sep 2026: the block appears only when the engine ran the check (c.web) ---- */
  if(c.web) h+=msbWebBlock(a,c,sec,nv.web); else nv.general.push(...nv.web);

  /* ---- 5.4 Compression Resistance N.b.Rd ---- */
  if(B && B.Fc>1e-9){
    h+=msbHead('Compression Resistance N.b.Rd');
    const NcrY=msbNcr(a.E,sec.Ix,B.LcrY), NcrZ=msbNcr(a.E,sec.Iy,B.LcrZ);
    const Ky=B.Ky, Kz=(B.KzEnd!=null? B.KzEnd : B.Kz);
    const kTag=B.leOverride? 'user L<sub>E</sub>' : (B.cantStrut? 'cantilever 2.0L' : 'end fixities [verify]');
    h+=msbRow('L<sub>ey</sub> = K<sub>y</sub>.L<sub>y</sub>', g(Ky,2)+' x '+msbM(S.L)+' = ('+msbEsc(B.leOverride? B.lcrBasis : B.lcrBasisY)+')', msbM(B.LcrY/1000)+' m', kTag);
    const Acm=(B.aeffOn? B.Aeff : sec.A*100)/100, Albl=B.aeffOn? 'A<sub>eff</sub>' : 'A';
    if(B.aeffOn) h+=msbRow('A<sub>eff</sub> (Class-4 web in uniform compression)', 'N<sub>Rk</sub> = A<sub>eff</sub>f<sub>y</sub>; &lambda;&#772; = &radic;(A<sub>eff</sub>f<sub>y</sub>/N<sub>cr</sub>) (6.3.1.3(1)); Table 6.7 Class-4 column with W<sub>eff.y</sub> = W<sub>el.y</sub> (flanges Class &le; 3; web Class 4 only in uniform compression, e<sub>N</sub> = 0) [assumption printed]', f1(Acm,2)+' cm&sup2;', 'EN 1993-1-5 4.4');
    h+=msbRow('&lambda;&#772;<sub>y</sub> = &radic;'+Albl+'.f<sub>y</sub>/N<sub>cr</sub>', '&radic;'+f1(Acm,2)+'x'+msbInt(c.fy)+'/'+f1(NcrY,2)+' (N<sub>cr,y</sub> = &pi;&sup2;EI<sub>y</sub>/L<sub>ey</sub>&sup2;)', msbR(B.lamY), '');
    h+=msbRow('N<sub>b.y.Rd</sub> = '+(B.aeffOn? 'A<sub>eff</sub>' : 'Area')+'.&chi;.f<sub>y</sub>/&gamma;<sub>M1</sub>', f1(Acm,2)+'x'+msbR(B.chiY)+'x'+msbInt(c.fy)+'/10/1 =', msbKN(B.NbY)+' kN', 'Curve '+B.cvY.curve);
    h+=msbRow('L<sub>ez</sub> = K<sub>z</sub>.L<sub>z</sub>', B.lczFromRestraints? 'largest lateral-restraint spacing (&le; '+g(Kz,2)+' x '+msbM(S.L)+' from the end fixities) =' : g(Kz,2)+' x '+msbM(S.L)+' = ('+msbEsc(B.leOverride? B.lcrBasis : B.lcrBasisZ)+')', msbM(B.LcrZ/1000)+' m', B.lczFromRestraints? 'P360 6.2' : kTag);
    h+=msbRow('&lambda;&#772;<sub>z</sub> = &radic;'+Albl+'.f<sub>y</sub>/N<sub>crz</sub>', '&radic;'+f1(Acm,2)+'x'+msbInt(c.fy)+'/'+f1(NcrZ,2)+' (N<sub>cr,z</sub> = &pi;&sup2;EI<sub>z</sub>/L<sub>ez</sub>&sup2;)', msbR(B.lamZ), '');
    h+=msbRow('N<sub>b.z.Rd</sub> = '+(B.aeffOn? 'A<sub>eff</sub>' : 'Area')+'.&chi;.f<sub>y</sub>/&gamma;<sub>M1</sub>', f1(Acm,2)+'x'+msbR(B.chiZ)+'x'+msbInt(c.fy)+'/10/1 =', msbKN(B.NbZ)+' kN', 'Curve '+B.cvZ.curve);
    // [beam-v03 addition, 19 Sep 2026 G3] channel: torsional / torsional-flexural buckling (cl 6.3.1.4)
    if(B.tfb && B.tfb.ok){
      const t=B.tfb;
      h+=msbRow('i<sub>0</sub>&sup2; = i<sub>y</sub>&sup2; + i<sub>z</sub>&sup2; + y<sub>0</sub>&sup2;', f1(t.iy,2)+'&sup2; + '+f1(t.iz,2)+'&sup2; + '+f1(t.y0,2)+'&sup2; (y<sub>0</sub> = '+msbEsc(t.y0Src)+')', f1(t.i0sq,1)+' mm&sup2;', '6.3.1.4');
      h+=msbRow('N<sub>cr.T</sub> = (GI<sub>T</sub> + &pi;&sup2;EI<sub>w</sub>/L<sub>T</sub>&sup2;)/i<sub>0</sub>&sup2;', '(81000 x '+g(t.IT/1e4,2)+'e4 + &pi;&sup2; x '+msbInt(a.E)+' x '+g(t.Iw/1e12,5)+'e12/'+g(t.LT,0)+'&sup2;)/'+f1(t.i0sq,1)+' (I<sub>T</sub>, I<sub>w</sub>: '+msbEsc(t.ITSrc)+'; L<sub>T</sub> = '+msbM(t.LT/1000)+' m = '+msbEsc(t.LTSrc)+')', msbKN(t.NcrT)+' kN', '6.3.1.4');
      h+=msbRow('N<sub>cr.TF</sub> = (N<sub>cr.y</sub> + N<sub>cr.T</sub>)/2&beta; [1 &minus; &radic;(1 &minus; 4&beta;N<sub>cr.y</sub>N<sub>cr.T</sub>/(N<sub>cr.y</sub> + N<sub>cr.T</sub>)&sup2;)]', 'N<sub>cr.y</sub> = '+msbKN(t.NcrY)+' (flexure about the axis of symmetry y-y, L<sub>ey</sub>), &beta; = 1 &minus; (y<sub>0</sub>/i<sub>0</sub>)&sup2; = '+msbR(t.beta)+'; N<sub>cr</sub> = min(N<sub>cr.T</sub>, N<sub>cr.TF</sub>) = '+msbKN(t.Ncr)+' ('+t.mode+' mode)', msbKN(t.NcrTF)+' kN', '6.3.1.4');
      h+=msbRow('&lambda;&#772;<sub>T</sub> = &radic;A.f<sub>y</sub>/N<sub>cr</sub>', '&radic;'+f1(sec.A,2)+'x'+msbInt(c.fy)+'/'+msbKN(t.Ncr), msbR(t.lamT), '');
      h+=msbRow('N<sub>b.T.Rd</sub> = Area.&chi;<sub>T</sub>.f<sub>y</sub>/&gamma;<sub>M1</sub>', f1(sec.A,2)+'x'+msbR(t.chiT)+'x'+msbInt(c.fy)+'/10/1 = (curve related to z-z: Table 6.2 U-sections, any axis [verify])', msbKN(t.NbT)+' kN', 'Curve '+t.cvT.curve);
      h+=msbRow('N<sub>Ed</sub>/N<sub>b.T.Rd</sub>', msbKN(B.Fc)+' / '+msbKN(t.NbT)+' = (the lower of &chi;<sub>T</sub> and the flexural &chi; feeds U<sub>N.y</sub>, U<sub>N.z</sub>)', msbR(t.util), msbWarn(t.util<=1.0001));
    } else {
      // 20 Sep 2026: reference software prints a torsional-buckling strut line (L_et, lambda_T, N_b.T.Rd) for every section; beam-v03 evaluates cl 6.3.1.4 for channels only
      h+=msbRow('L<sub>et</sub> = K<sub>t</sub>.L<sub>z</sub> ; &lambda;&#772;<sub>T</sub> ; N<sub>b.T.Rd</sub>', 'n/a - not evaluated by beam-v03 (cl 6.3.1.4 torsional / torsional-flexural buckling is evaluated for channels only)', '&mdash;', 'not evaluated');
    }
    {
      // N.Ed/N.b.Rd against the lower flexural resistance (the values U_N.y / U_N.z use; N_b.T.Rd already folded in for a channel)
      const NbYe=(B.NbYeff!=null? B.NbYeff : B.NbY), NbZe=(B.NbZeff!=null? B.NbZeff : B.NbZ), NbMin=Math.min(NbYe,NbZe), uN=B.Fc/Math.max(NbMin,1e-9);
      h+=msbRow('N<sub>Ed</sub>/N<sub>b.Rd</sub>', msbKN(B.Fc)+' / '+msbKN(NbMin)+' = (min of N<sub>b.y.Rd</sub>, N<sub>b.z.Rd</sub>'+((B.tfb&&B.tfb.ok)? ', N<sub>b.T.Rd</sub>' : '')+')', msbR(uN), msbWarn(uN<=1.0001));
    }
    h+=msbNotVerifiedRows(nv.compression);
  }

  /* ---- 5.5 Equivalent Uniform Moment Factor(s) ---- */
  // 20 Sep 2026 review: one rule on both Mcr routes - reference software prints "C1 = ... Uniform" for a hollow section too
  // (SHS-L2 / RHS-L2 printouts), so the standard-route box with lambda_LT <= 0.4 (LT.na) no longer suppresses the block
  const c1Line=!fullRest && LT && !LT.failed;
  if(c1Line || B){
    // 20 Sep 2026: the C_mLT / C_mz / C_my lines exist only when the cl 6.3.3 interaction is evaluated (B), so the heading follows B, not the brief type
    h+=msbHead(B? 'Equivalent Uniform Moment Factors C1, C.mLT, C.mz, and C.my' : 'Equivalent Uniform Moment Factor C1');
    if(c1Line){
      const ci=LT.c1in;
      const ciTxt= ci? f1(ci.M1,1)+', '+f1(ci.M2,1)+', '+f1(ci.Mo,1)+', '+f1(ci.psi,3)+', '+f1(ci.mu,3) : '&mdash;';
      if(eigen){
        const lbl=msbEsc(LT.c1label);
        if(S.C1o!=null) h+=msbRow('C<sub>1</sub> = user override', 'user override; '+lbl, msbR(LT.C1), 'User');
        else if(/critical bay/.test(lbl)) h+=msbRow('C<sub>1</sub> = M<sub>cr</sub>/M<sub>cr,uniform</sub> (critical bay)', 'bay between restraints containing the eigenmode peak (x = '+msbM((LT.modePeakX||0)/1000)+' m); whole-member ratio '+msbR(LT.McrShape/Math.max(LT.McrUniform,1e-9)), msbR(LT.C1), 'NA 2.18');
        else h+=msbRow('C<sub>1</sub> = M<sub>cr</sub>/M<sub>cr,uniform</sub>', msbKNm(LT.McrShape)+' / '+msbKNm(LT.McrUniform)+' kN.m (shape only: z<sub>g</sub> = z<sub>j</sub> = 0)', msbR(LT.C1), LT.c1Trusted===false? 'not converged: k<sub>c</sub> = 1' : 'k<sub>c</sub> only');
        if(LT.std && LT.std.C1!=null && isFinite(LT.std.C1)){
          const sg=LT.c1seg||{};
          const segTxt=(sg.whole===false)? ' [segment '+msbM(sg.xa/1000)+'&ndash;'+msbM(sg.xb/1000)+' m]' : '';
          h+=msbRow((LT.std.route==='sn006a'? 'C' : 'C<sub>1</sub>')+' = fn(M<sub>1</sub>, M<sub>2</sub>, M<sub>o</sub>, &psi;, &mu;) (standard, comparison)', ciTxt+segTxt+' &mdash; '+msbEsc(LT.std.label), msbR(LT.std.C1), msbC1Tag(LT.std.label,LT.std.route));
        }
      } else {
        const C1v=(LT.C1show!=null? LT.C1show : c.C1);
        const tag= LT.cant? 'Cantilever' : msbC1Tag(LT.c1label||c.c1label, LT.c1route);
        // 20 Sep 2026: the cantilever prints the reference software SN006a form "C1 = fn(M, Zg, kwt) ... Ncci-sn006" (LT.kwt, LT.zg, LT.warp of the engine)
        // 20 Sep 2026 review: a channel on the P362 kappa chain uses no C1 (lambda_LT = (L_e/i_z)/kappa) unless the
        // shear-centre Mcr route (LT.chanMcr, SN003a form) is offered; a channel cantilever / torsion case prints the
        // table value as "not used" rather than an I-section "Cantilever [SN003a]" tag
        const chanNoC1=!!(LT.channel && !LT.chanMcr);
        if(LT.cant) h+=msbRow('C<sub>1</sub> = fn(M, Z<sub>g</sub>, &kappa;<sub>wt</sub>)', 'M<sub>1</sub>, M<sub>2</sub> = '+(ci? f1(ci.M1,1)+', '+f1(ci.M2,1) : '&mdash;')+' kN.m; z<sub>g</sub> = '+g(LT.zg||0,0)+' mm; &kappa;<sub>wt</sub> = '+msbDash(LT.kwt,v=>f1(v,3))+'; "Cantilever end warping '+(LT.warp==='restr'? 'fixed' : 'free')+'"', msbR(C1v), 'Ncci-sn006');
        else if(chanNoC1) h+=msbRow('C<sub>1</sub> = fn(M<sub>1</sub>, M<sub>2</sub>, M<sub>o</sub>, &psi;, &mu;)', ciTxt+' &mdash; not used: the P362 &kappa; chain sets &lambda;&#772;<sub>LT</sub> = (L<sub>e</sub>/i<sub>z</sub>)/&kappa; without C<sub>1</sub>'+(isCant? ' (channel cantilever: the SN006a coefficients are for I sections)' : ''), '&mdash;', 'not used');
        else h+=msbRow('C<sub>1</sub> = fn(M<sub>1</sub>, M<sub>2</sub>, M<sub>o</sub>, &psi;, &mu;)', ciTxt, msbR(C1v), tag);
        h+=msbRow('C<sub>1</sub> basis', msbEsc(LT.c1label||c.c1label), '', LT.cant? 'SN006a' : chanNoC1? 'P362' : 'SN003a', 'ms-basis');
      }
    }
    if(B){
      const aB=(B.combo && nULS>1)? a.ulsResults.find(r=>r.combo.label===B.combo) : null;
      const aCm= aB? analysisForCombination(a,aB) : a;
      const p=msbCmB3Params(aCm);
      const pTxt='M<sub>h</sub> = '+f1(p.Mh,2)+', M<sub>s</sub> = '+f1(p.Ms,2)+', &psi; = '+f1(p.psi,3)+', '+(p.alphaS!=null? '&alpha;<sub>s</sub> = '+f1(p.alphaS,3) : '&alpha;<sub>h</sub> = '+f1(p.alphaH,3));
      const form=msbCmB3Form(B.cmLabel,p.alphaS)+(B.swayNote? ' &ge; 0.9 (sway mode)' : '');
      h+=msbRow('C<sub>mLT</sub> = '+form, pTxt, msbR(B.CmLT), 'Table B.3');
      // 20 Sep 2026 torsion + N/Mz: with torsion M_z,Ed of the cl 6.3.3 interaction = imposed M_z + max |phi.M_y| (annexB2)
      // 20 Sep 2026 review: psi = 1 describes the imposed constant M_z only; the twist part phi(x)M_y(x) is a hump diagram taken at the Table B.3 upper bound C_mz = 1.0
      h+=msbRow('C<sub>mz</sub> = Max(0.6+0.4&psi;, 0.4)', (T && B.MzTwist>1e-9? 'M<sub>z.Ed</sub> = '+msbKNm(B.MzImp)+' imposed + '+msbKNm(B.MzTwist)+' twist &phi;.M<sub>y</sub> = '+msbKNm(B.MzEd)+'; &psi; = 1 for the imposed constant M<sub>z</sub>, the twist part taken at C<sub>mz</sub> = 1.0 (Table B.3 upper bound, conservative)' : 'M = '+msbKNm(B.MzEd)+', &psi; = 1.000'), msbR(B.Cmz), 'Table B.3');
      h+=msbRow('C<sub>my</sub> = '+form, pTxt, msbR(B.Cmy), 'Table B.3');
    }
  }

  /* ---- 5.6 Lateral Buckling Check M.b.Rd ---- */
  h+=msbHead('Lateral Buckling Check M.b.Rd');
  const WyTxt=f1(Wy_cm3,1);
  const lamLine=(lam,Mcr)=>msbRow('&lambda;&#772;<sub>LT</sub> = &radic;W.f<sub>y</sub>/M<sub>cr</sub>', '&radic; '+WyTxt+' x '+msbInt(c.fy)+' / '+msbKNm(Mcr), msbR(lam), '');
  const ignLine=(lam)=>msbRow('&lambda;&#772;<sub>LT</sub> &le; &lambda;&#772;<sub>LT,0</sub>', msbR(lam)+' &le; 0.4 : LTB may be ignored', '&chi;<sub>LT</sub> = 1.000', '6.3.2.2(4)');
  const chiLine=(lam,Phi,alphaLT,chi,curve)=>msbRow('&chi;<sub>LT</sub> = Fn(&lambda;&#772;<sub>LT</sub>, &Phi;<sub>LT</sub>, &alpha;<sub>LT</sub>)', msbR(lam)+', '+msbDash(Phi,v=>f1(v,3))+', '+g(alphaLT,2), msbR(chi), 'Curve '+curve);
  const chiModLine=(chi,lam,kc,f,chiMod,tag)=>msbRow('&chi;<sub>LT.mod</sub> = Fn(&chi;<sub>LT</sub>, &lambda;&#772;<sub>LT</sub>, k<sub>c</sub>, f)', msbR(chi)+', '+msbR(lam)+', '+msbDash(kc,v=>f1(v,3))+', '+msbR(f), msbR(chiMod), tag);
  const mbLine=(chiMod,Mb,capped)=>msbRow('M<sub>b.Rd</sub> = &chi;'+Wlbl+'.f<sub>y</sub>'+(capped!==false? ' &le; M<sub>c.y.Rd</sub>' : ''), msbR(chiMod)+' x '+WyTxt+' x '+msbInt(c.fy)+(capped!==false? ' &le; '+msbKNm(c.McRd) : '')+' =', msbKNm(Mb)+' kN.m', '');
  const ratioLine=(Mx,Mb)=>msbRow('M<sub>y.Ed</sub>/M<sub>b.Rd</sub>', msbKNm(Mx)+' / '+msbKNm(Mb), (Mb>0&&isFinite(c.ltbUtil))? msbR(c.ltbUtil) : '&mdash;', (Mb>0&&isFinite(c.ltbUtil))? msbWarn(c.ltbUtil<=1.0001) : 'not evaluated');
  // 20 Sep 2026 review: the reference software "Section not susceptible to lateral torsional buckling" line of a hollow section
  // with lambda_LT <= 0.4, the same on both Mcr routes (the eigen route printed the chi_LT = 1 / chi_LT.mod pair instead)
  const boxIgnRow=(lam,Mb)=> S.family==='shs'
    ? msbRow('M<sub>b.Rd</sub> = M<sub>c.y.Rd</sub>', 'closed hollow section &mdash; not susceptible to LTB (&lambda;&#772;<sub>LT</sub> = '+msbR(lam)+' &le; 0.4)', msbKNm(Mb)+' kN.m', '6.3.2.1(2)')
    : msbRow('M<sub>b.Rd</sub> = M<sub>c.y.Rd</sub>', 'closed hollow section, &lambda;&#772;<sub>LT</sub> = '+msbR(lam)+' &le; 0.4: LTB may be ignored', msbKNm(Mb)+' kN.m', '6.3.2.2(4)');
  const leK=ltbLeFactor()*(S.destab?1.2:1);
  if(fullRest){
    h+=msbRow('M<sub>b.Rd</sub> = M<sub>c.y.Rd</sub>', 'Fully Restrained', msbKNm(c.McRd)+' kN.m', '');
  } else if(!LT){
    h+=msbRow('M<sub>b.Rd</sub>', 'no LTB result in the check object', '&mdash;', 'not evaluated');
  } else if(LT.box && !eigen){
    // closed section on the standard route: SN003a chain with I_w = 0
    if(LT.ignM){
      h+=boxIgnRow(LT.lamLTmcr,LT.MbRd);
      if(LT.Mcr>0) h+=msbRow('M<sub>cr</sub> (information, I<sub>w</sub> = 0)', 'SN003a with C<sub>1</sub> = '+msbR(c.C1)+', L<sub>e</sub> = '+msbM(c.LE/1000)+' m; '+msbEsc(LT.zgNote||''), msbKNm(LT.Mcr)+' kN.m', 'SN003a');
    } else {
      h+=msbRow('L<sub>e</sub> = '+g(leK,2)+' L', g(leK,2)+' x '+msbM(S.L)+' =', msbM(c.LE/1000)+' m', '');
      h+=msbRow('M<sub>cr</sub> = Fn(C<sub>1</sub>, L<sub>e</sub>, I<sub>z</sub>, I<sub>t</sub>, I<sub>w</sub> = 0, E)', msbR(c.C1)+', '+f1(c.LE/1000,3)+', '+g(sec.Iy,2)+', '+g(sec.J,2)+', 0, '+msbInt(a.E)+'; '+msbEsc(LT.zgNote||''), msbKNm(LT.Mcr)+' kN.m', LT.zgBlocked? '<span class="ms-warn">BLOCKED</span>' : 'SN003a');
      h+=lamLine(LT.lamLTmcr,LT.Mcr);
      h+=chiLine(LT.lamLTmcr,LT.PhiM,LT.curve.alphaLT,LT.chiM,LT.curve.curve);
      h+=chiModLine(LT.chiM,LT.lamLTmcr,LT.kc,LT.fM,LT.chiModM,'6.3.2.3 / NA Table NA.1');
      h+=mbLine(LT.chiModM,LT.MbRd);
    }
    h+=ratioLine(c.Mx,LT.MbRd);
  } else if(eigen){
    if(LT.failed){
      h+=msbRow('M<sub>cr</sub> = FE eigenvalue', msbEsc(LT.err), '&mdash;', '<span class="ms-warn">BLOCKED</span>');
      h+=ratioLine(c.Mx,0);
    } else {
      const sg=LT.spanGoverns? LT.spanGov : null;
      const endBc=ltbEndBcText(S, sec, true);
      if(sg) h+=msbRow('Governing bay '+msbM(sg.a/1000)+'&ndash;'+msbM(sg.b/1000)+' m (isolated, fork ends)', 'M<sub>y.Ed</sub> = '+msbKNm(sg.Ms)+', M<sub>cr</sub> = '+msbKNm(sg.Mcr)+', &lambda;&#772; = '+msbR(sg.lam)+', &chi; = '+msbR(sg.chi), msbKNm(sg.Mb)+' kN.m', 'bay by bay');
      h+=msbRow('L<sub>e</sub> = portion between restraints', endBc+((LT.vPoints||[]).length>2? '; lateral restraint points at x = '+(LT.vPoints||[]).map(x=>msbM(x/1000)).join(', ')+' m' : ''), msbM((xb-xa)/1000)+' m', 'FE');
      const zgTxt=(LT.zgValues&&LT.zgValues.length>1)? '; z<sub>g</sub> = '+LT.zgValues.map(z=>g(z,0)).join(', ')+' mm' : (Math.abs(LT.zg||0)>1e-9? '; z<sub>g</sub> = '+g(LT.zg,0)+' mm (load reversed: '+msbKNm(LT.McrRev)+')' : '');
      const mcrBad = LT.mcrConverged===false || (LT.meshError||0)>0.005;
      h+=msbRow('M<sub>cr</sub> = FE eigenvalue (n<sub>Elem</sub>, mesh error)', msbInt(LT.nElem)+' elements, '+f1((LT.meshError||0)*100,3)+' %'+(LT.nCombos>1? '; governing: '+msbEsc(LT.governCombo) : '')+zgTxt+(LT.nSolves!=null? '; '+LT.nCombos+' combination(s), '+LT.nSolves+' solve(s)'+(LT.nCached? ' + '+LT.nCached+' cached' : '') : ''), msbKNm(LT.Mcr)+' kN.m'+(sg? ' (whole member)' : ''), mcrBad? '<span class="ms-warn">BLOCKED</span>' : 'converged');
      const lam= sg? sg.lam : LT.lamLT, Mcr= sg? sg.Mcr : LT.Mcr;
      h+=lamLine(lam,Mcr);
      if(sg){
        if(lam<=0.4) h+=ignLine(lam);
        else h+=chiLine(lam,msbPhiLT(lam,LT.curve.alphaLT),LT.curve.alphaLT,sg.chi,LT.curve.curve);
        h+=chiModLine(sg.chi,lam,null,1,sg.chi,'isolated bay: f = 1');
        h+=mbLine(sg.chi,sg.Mb);
        h+=ratioLine(sg.Ms,sg.Mb);
      } else if(LT.box && (LT.ign || lam<=0.4)){
        h+=boxIgnRow(lam,LT.MbRd);   // the FE Mcr and lambda rows above stay as information
        h+=ratioLine(LT.MxGov!=null? LT.MxGov : c.Mx, LT.MbRd);
      } else {
        if(LT.ign) h+=ignLine(lam);
        else h+=chiLine(lam,LT.Phi,LT.curve.alphaLT,LT.chi,LT.curve.curve);
        h+=chiModLine(LT.chi,lam,LT.kc,LT.f,LT.chiMod,isCant? 'f = 1 (cantilever)' : (LT.kcFloored? '6.3.2.3; k<sub>c</sub> floored at 0.60 (Table 6.6)' : '6.3.2.3'));
        h+=mbLine(LT.chiMod,LT.MbRd);
        h+=ratioLine(LT.MxGov!=null? LT.MxGov : c.Mx, LT.MbRd);
      }
    }
  } else if(LT.cant){
    h+=msbRow('M<sub>cr0</sub> = (&pi;/L)&radic;(EI<sub>z</sub>GI<sub>t</sub>)', 'L = '+msbM(S.L)+' m, E = '+msbInt(a.E)+', I<sub>z</sub> = '+g(sec.Iy,2)+' cm&#8308;, I<sub>t</sub> = '+g(sec.J,2)+' cm&#8308;, G = 81000', msbKNm(LT.Mcr0)+' kN.m', 'SN006a');
    const blocked = !(LT.C>0);
    h+=msbRow('C = Fn(&kappa;<sub>wt</sub>, &eta;, warping)', f1(LT.kwt,3)+', '+f1(LT.eta,3)+', '+(LT.warp==='restr'? 'restrained' : 'free')+'; '+msbEsc(LT.caseLbl)+(LT.Cq!=null&&LT.CF!=null? '; C<sub>q</sub> = '+f1(LT.Cq,3)+', C<sub>F</sub> = '+f1(LT.CF,3)+', Eq (7)' : ''), blocked? '&mdash;' : msbR(LT.C), blocked? '<span class="ms-warn">BLOCKED</span>' : 'Tables 3.1&ndash;3.3');
    if(!blocked){
      h+=msbRow('M<sub>cr</sub> = C&middot;M<sub>cr0</sub>', msbR(LT.C)+' x '+msbKNm(LT.Mcr0), msbKNm(LT.Mcr)+' kN.m', '');
      h+=lamLine(LT.lamLTmcr,LT.Mcr);
      if(LT.ignM) h+=ignLine(LT.lamLTmcr); else h+=chiLine(LT.lamLTmcr,LT.PhiM,LT.curve.alphaLT,LT.chiM,LT.curve.curve);
      h+=chiModLine(LT.chiM,LT.lamLTmcr,null,1,LT.chiModM,'f = 1 (cantilever)');
      h+=mbLine(LT.chiModM,LT.MbRd);
    }
    h+=ratioLine(c.Mx,LT.MbRd);
  } else if(LT.channel){
    h+=msbRow('L<sub>e</sub> = '+g(leK,2)+' L', g(leK,2)+' x '+msbM(S.L)+' =', msbM(c.LE/1000)+' m', '');
    if(Math.abs(LT.zg||0)>1e-9) h+=msbRow('z<sub>g</sub>', msbEsc(LT.zgNote||''), '', 'load height');
    h+=msbRow('&lambda;&#772;<sub>LT</sub> = (L<sub>e</sub>/i<sub>z</sub>)/&kappa;', '('+g(c.LE,0)+'/'+g(LT.ry,1)+')/'+g(LT.kappa,0)+' ('+S.grade+')', msbR(LT.lamLTmcr), 'P362 channel');
    if(LT.ignM) h+=ignLine(LT.lamLTmcr); else h+=chiLine(LT.lamLTmcr,LT.PhiM,0.76,LT.chiM,'d');
    h+=msbRow('M<sub>b.Rd</sub> = &chi;W<sub>pl.y</sub>.f<sub>y</sub> (&kappa; chain, no f)', msbR(LT.chiM)+' x '+WyTxt+' x '+msbInt(c.fy)+' =', msbKNm(LT.MbSimp)+' kN.m', '');
    h+=msbRow('M<sub>cr</sub> (back-calculated) = W.f<sub>y</sub>/&lambda;&#772;<sub>LT</sub>&sup2;', WyTxt+' x '+msbInt(c.fy)+' / '+msbR(LT.lamLTmcr)+'&sup2;', msbKNm(LT.McrBack)+' kN.m', 'Annex A only');
    if(LT.chanMcr){
      const cm=LT.chanMcr;
      h+=msbRow('M<sub>cr</sub> route (load through the shear centre)', 'M<sub>cr</sub> = '+msbKNm(cm.Mcr)+', &lambda;&#772; = '+msbR(cm.lam)+', &chi; = '+msbR(cm.chi)+', f = '+msbR(cm.f)+', &chi;<sub>mod</sub> = '+msbR(cm.chiMod), msbKNm(cm.Mb)+' kN.m', 'SN003a (z<sub>j</sub> = 0)');
    }
    if(Math.abs(LT.MbRd-LT.MbSimp)>1e-6) h+=msbRow('M<sub>b.Rd</sub> (design basis)', msbEsc(c.ltbBasis), msbKNm(LT.MbRd)+' kN.m', '');
    h+=ratioLine(c.Mx,LT.MbRd);
  } else {
    // D. standard closed form, SN003a (I/H)
    h+=msbRow('L<sub>e</sub> = '+g(leK,2)+' L', g(leK,2)+' x '+msbM(S.L)+' =', msbM(c.LE/1000)+' m', '');
    h+=msbRow('M<sub>cr</sub> = Fn(C<sub>1</sub>, L<sub>e</sub>, I<sub>z</sub>, I<sub>t</sub>, I<sub>w</sub>, E)', msbR(c.C1)+', '+f1(c.LE/1000,3)+', '+g(sec.Iy,2)+', '+g(sec.J,2)+', '+g(sec.Iw||0,4)+', '+msbInt(a.E)+(LT.zgUsed? ', C<sub>2</sub>z<sub>g</sub> = '+g(LT.C2*LT.zg,1)+' mm' : '')+'; '+msbEsc(LT.zgNote||''), msbKNm(LT.Mcr)+' kN.m', LT.zgBlocked? '<span class="ms-warn">BLOCKED</span>' : 'SN003a');
    h+=lamLine(LT.lamLTmcr,LT.Mcr);
    if(LT.ignM) h+=ignLine(LT.lamLTmcr); else h+=chiLine(LT.lamLTmcr,LT.PhiM,LT.curve.alphaLT,LT.chiM,LT.curve.curve);
    h+=chiModLine(LT.chiM,LT.lamLTmcr,LT.kc,LT.fM,LT.chiModM,LT.kcFloored? '6.3.2.3; k<sub>c</sub> floored at 0.60 (Table 6.6)' : '6.3.2.3');
    h+=mbLine(LT.chiModM,LT.MbRd);
    // [beam-v03 addition] the P362 Expn 6.55 simplified slenderness, comparison only (never the design basis)
    h+=msbRow('&lambda;&#772;<sub>LT</sub> (P362 6.55 simplified) ; M<sub>b.Rd</sub> (comparison only)', '(1/&radic;C<sub>1</sub>)0.9&lambda;&#772;<sub>z</sub>&radic;&beta;<sub>w</sub> = '+msbR(LT.lamLTsimp)+'; &chi;<sub>LT.mod</sub> = '+msbR(LT.chiModS)+'; M<sub>y.Ed</sub>/M<sub>b.Rd,simp</sub> = '+msbR(c.Mx/Math.max(LT.MbSimp,1e-9)), msbKNm(LT.MbSimp)+' kN.m', 'P362 6.55 comparison');
    h+=ratioLine(c.Mx,LT.MbRd);
  }
  h+=msbNotVerifiedRows(nv.ltb);
  // 20 Sep 2026 review: the eigen route's LTB notes (LT.warn: L_E factor ignored on this route, warping flag not
  // applied to a closed section, C1 not trusted / k_c floored, L_cr,z from the restraint spacing) were printed only
  // by the deleted CED report; they are advisory rows of this block now
  h+=msbAdvisoryRows(LT && LT.warn);
  /* ---- 7.2 portion table ---- */
  const RF=c.restraintForces||null;
  const rfRows=(RF && RF.rows)||[];
  const hasSeg=!!(LT && LT.segments && LT.segments.length);
  // 20 Sep 2026 review: a "Lateral Restraint Portions" block only with intermediate lateral restraints (bay by bay
  // and / or the restraint-force table), as reference software prints no portion for a plain cantilever or fixed-ended
  // beam; the restraint design force of a support restraint with M_Ed != 0 (fixed end, cantilever root) is an
  // advisory row of this block instead. A simply supported beam (M_Ed = 0 at both supports) prints nothing.
  const rfInter=rfRows.some(r=> r.kind==='lateral' || (r.x>1e-6 && r.x<a.L-1e-6));
  if(RF && !rfInter && !hasSeg){
    rfRows.filter(r=>Math.abs(r.F)>5e-4).forEach(r=>{   // 5e-4 kN: below the printed 3 dp (a pinned end carries ~1e-5 kN.m numerically)
      h+=msbRow('Restraint design force @ '+msbM(r.x/1000)+' m', msbEsc(r.label)+': M<sub>Ed</sub> = '+msbKNm(r.MEd)+' kN.m ('+msbEsc(r.combo)+'); N<sub>f.Ed</sub> = M<sub>Ed</sub>/h = '+msbKN(r.NfEd)+' kN (h = '+msbMM(RF.h)+' mm); 2.5 % N<sub>f.Ed</sub> = (restraint design force, advisory: 6.3.5.2(5)(b), not part of the verdict)', msbKN(r.F)+' kN', 'advisory', 'ms-advrow');
    });
  }
  h+=msbModeShape(LT);   // 20 Sep 2026 review: the eigen route's buckled mode shape (was in the deleted CED report)
  if(hasSeg || rfInter){
    const rfShow=rfInter;
    h+=msbHead(hasSeg? 'Lateral Restraint Portions (bay by bay, fork ends)' : 'Lateral Restraint Portions (restraint design forces)');
    if(hasSeg){
      let worst=null; LT.segments.forEach(s2=>{ if(s2.ok && (!worst||s2.util>worst.util)) worst=s2; });
      h+='<table class="ms-combos"><thead><tr><th>Portion</th><th>From &ndash; To (m)</th><th>L<sub>e</sub> (m)</th><th>M<sub>y.Ed</sub> (kN.m)</th><th>M<sub>cr</sub> (kN.m)</th><th>&lambda;&#772;<sub>LT</sub></th><th>&chi;<sub>LT</sub></th><th>M<sub>b.Rd</sub> (kN.m)</th><th>M<sub>y.Ed</sub>/M<sub>b.Rd</sub></th><th></th></tr></thead><tbody>'+
         LT.segments.map((s2,i)=>'<tr><td class="num">'+(i+1)+'</td><td class="num">'+msbM(s2.a/1000)+' &ndash; '+msbM(s2.b/1000)+'</td><td class="num">'+msbM((s2.b-s2.a)/1000)+'</td>'+
           (s2.ok? '<td class="num">'+msbKNm(s2.Ms)+'</td><td class="num">'+msbKNm(s2.Mcr)+'</td><td class="num">'+msbR(s2.lam)+'</td><td class="num">'+msbR(s2.chi)+'</td><td class="num">'+msbKNm(s2.Mb)+'</td><td class="num">'+msbR(s2.util)+'</td><td>'+(s2===worst? (LT.spanGoverns? 'governs' : 'worst segment') : '')+'</td>'
                  : '<td colspan="7">not solved: '+msbEsc(s2.err)+'</td>')+'</tr>').join('')+
         '</tbody></table><div class="ms-note">&chi;<sub>LT</sub> per portion without the f-factor. Whole-member M<sub>cr</sub> = '+msbKNm(LT.Mcr)+' kN.m. '+msbEsc(c.ltbBasis)+'</div>';
    }
    // [beam-v03 addition, 19 Sep 2026 G3] restraint design forces (advisory): 2.5 % of N_f,Ed = M_Ed/h at every restraint station
    if(rfShow){
      h+='<table class="ms-combos ms-restraint"><thead><tr><th>x (m)</th><th>Restraint</th><th>M<sub>Ed</sub> (kN.m)</th><th>Load case</th><th>N<sub>f.Ed</sub> = M<sub>Ed</sub>/h (kN)</th><th>2.5 % N<sub>f.Ed</sub> (kN)</th><th></th></tr></thead><tbody>'+
        RF.rows.map(r=>'<tr><td class="num">'+msbM(r.x/1000)+'</td><td>'+msbEsc(r.label)+'</td><td class="num">'+msbKNm(r.MEd)+'</td><td>'+msbEsc(r.combo)+'</td><td class="num">'+msbKN(r.NfEd)+'</td><td class="num">'+msbKN(r.F)+'</td><td>restraint design force, advisory</td></tr>').join('')+
        '</tbody></table><div class="ms-note">h = '+msbMM(RF.h)+' mm. '+RF.basis+'</div>';
    }
  }

  /* ---- 5.7 Buckling Resistance ---- */
  // 20 Sep 2026: reference software prints this block in every Axial with Moments brief (zeros for N);
  // with a tensile N_Ed the engine evaluates no cl 6.3.3 interaction (B = null) and the block says so
  if(AX && !B){
    h+=msbHead('Buckling Resistance');
    h+=msbRow('U<sub>N.y</sub>, U<sub>N.z</sub>, U<sub>M.y</sub>, U<sub>M.z</sub>, k<sub>ij</sub>, Eq 6.61 / 6.62', 'n/a - not evaluated by beam-v03 (N<sub>Ed</sub> is tensile: cl 6.3.3 member buckling is not required; lateral-torsional buckling is verified above)', '&mdash;', 'not evaluated');
  }
  if(B && (B.Fc>1e-9 || B.biax)){
    h+=msbHead('Buckling Resistance');
    const MbEff=(B.MbRdEff!=null? B.MbRdEff : B.MbRdI);
    const UMy=B.Mx/Math.max(MbEff,1e-9);
    const nTag=(u)=> B.Fc>1e-9? msbWarn(u<=1.0001) : 'N<sub>Ed</sub> = 0';
    const NRk= B.aeffOn? B.Aeff*c.fy/1000 : msbNRk(sec,c.fy);
    const tfT=(B.tfb&&B.tfb.ok)? '; min with N<sub>b.T.Rd</sub> = '+msbKN(B.tfb.NbT)+' (6.3.1.4)' : '';
    if(B.aeffOn) h+=msbRow('Table 6.7 Class-4 column', 'N<sub>Rk</sub> = A<sub>eff</sub>f<sub>y</sub> = '+msbKN(NRk)+' kN; M<sub>y.Rk</sub> = W<sub>eff.y</sub>f<sub>y</sub> with W<sub>eff.y</sub> = W<sub>el.y</sub> = '+f1(sec.Zx,1)+' cm&sup3; (flanges Class &le; 3, web Class 4 only in uniform compression, &Delta;M = 0) [assumption]; k<sub>ij</sub> from the Class 3/4 rows', '', '6.3.3');
    h+=msbRow('U<sub>N.y</sub> = N<sub>Ed</sub>/(&chi;<sub>y</sub>.N<sub>Rk</sub>/&gamma;<sub>M1</sub>)', msbKN(B.Fc)+' / '+msbKN(B.NbYeff!=null? B.NbYeff : B.NbY)+' (N<sub>Rk</sub> = '+msbKN(NRk)+tfT+')', msbR(B.ny), nTag(B.ny));
    h+=msbRow('U<sub>N.z</sub> = N<sub>Ed</sub>/(&chi;<sub>z</sub>.N<sub>Rk</sub>/&gamma;<sub>M1</sub>)', msbKN(B.Fc)+' / '+msbKN(B.NbZeff!=null? B.NbZeff : B.NbZ)+(tfT? ' ('+tfT.slice(2)+')' : ''), msbR(B.nz), nTag(B.nz));
    h+=msbRow('U<sub>M.y</sub> = M<sub>y.Ed</sub>/(&chi;<sub>LT</sub>.M<sub>y.Rk</sub>/&gamma;<sub>M1</sub>)', msbKNm(B.Mx)+' / '+msbKNm(MbEff)+' (M<sub>y.Rk</sub> = '+msbKNm(B.aeffOn? sec.Zx*1e3*c.fy/1e6 : msbMyRk(c.Wy,c.fy))+(B.aeffOn&&B.wFac!==1? '; W<sub>el.y</sub>/W<sub>pl.y</sub> = '+msbR(B.wFac)+' applied to M<sub>b.Rd</sub>' : '')+')', msbR(UMy), msbWarn(UMy<=1.0001));
    // 20 Sep 2026 torsion + N/Mz: M_z,Ed = imposed M_z + max |phi.M_y| of the torsion solution (annexB2, EN 1993-1-1 5.2.1(3)); printed as its two parts with torsion
    h+=msbRow('U<sub>M.z</sub> = M<sub>z.Ed</sub>/(M<sub>z.Rk</sub>/&gamma;<sub>M1</sub>)', (T && B.MzTwist>1e-9? '('+msbKNm(B.MzImp)+' imposed + '+msbKNm(B.MzTwist)+' twist &phi;.M<sub>y</sub> = '+msbKNm(B.MzEd)+')' : msbKNm(B.MzEd))+' / '+msbKNm(B.Mcz), msbR(B.mzTerm), msbWarn(B.mzTerm<=1.0001));
    // 20 Sep 2026: the reference software "kzy method" line (the engine's Table B.1 / B.2 choice)
    h+=msbRow('k<sub>zy</sub> method', B.useB1? 'not susceptible to torsional deformation (closed section, fully restrained, or M<sub>b.Rd</sub> = M<sub>c.y.Rd</sub>), using Table B.1' : 'M<sub>b.Rd</sub> &lt; M<sub>c.y.Rd</sub> therefore susceptible to LTB, using Table B.2', '', B.useB1? 'Table B.1' : 'Table B.2', 'ms-basis');
    if(B.c12){
      h+=msbRow('k<sub>yy</sub> = C<sub>my</sub>{1+(&lambda;&#772;<sub>y</sub>&minus;0.2)U<sub>N.y</sub>}', msbR(B.Cmy)+'{1+('+msbR(B.lamY)+'&minus;0.2)x'+msbR(B.ny)+'} &le; '+msbR(B.Cmy)+'(1+0.8x'+msbR(B.ny)+')', msbR(B.kyy), 'Table B.1');
      if(B.rhsRow){
        h+=msbRow('k<sub>zz</sub> = C<sub>mz</sub>{1+(&lambda;&#772;<sub>z</sub>&minus;0.2)U<sub>N.z</sub>}', msbR(B.Cmz)+'{1+('+msbR(B.lamZ)+'&minus;0.2)x'+msbR(B.nz)+'} &le; '+msbR(B.Cmz)+'(1+0.8x'+msbR(B.nz)+')', msbR(B.kzz), 'Table B.1 (RHS)');
        h+=msbRow('k<sub>yz</sub> = k<sub>zz</sub>', 'rectangular hollow section', msbR(B.kyz), 'Table B.1 (RHS)');
      } else {
        h+=msbRow('k<sub>zz</sub> = C<sub>mz</sub>{1+(2&lambda;&#772;<sub>z</sub>&minus;0.6)U<sub>N.z</sub>}', msbR(B.Cmz)+'{1+(2x'+msbR(B.lamZ)+'&minus;0.6)x'+msbR(B.nz)+'} &le; '+msbR(B.Cmz)+'(1+1.4x'+msbR(B.nz)+')', msbR(B.kzz), '');
        h+=msbRow('k<sub>yz</sub> = 0.6k<sub>zz</sub>', '0.6 x '+msbR(B.kzz), msbR(B.kyz), '');
      }
    } else {
      h+=msbRow('k<sub>yy</sub> = C<sub>my</sub>{1+0.6&lambda;&#772;<sub>y</sub>U<sub>N.y</sub>}', msbR(B.Cmy)+'{1+0.6x'+msbR(B.lamY)+'x'+msbR(B.ny)+'} &le; '+msbR(B.Cmy)+'(1+0.6x'+msbR(B.ny)+')', msbR(B.kyy), 'Table B.1 (Class 3)');
      h+=msbRow('k<sub>zz</sub> = C<sub>mz</sub>{1+0.6&lambda;&#772;<sub>z</sub>U<sub>N.z</sub>}', msbR(B.Cmz)+'{1+0.6x'+msbR(B.lamZ)+'x'+msbR(B.nz)+'} &le; '+msbR(B.Cmz)+'(1+0.6x'+msbR(B.nz)+')', msbR(B.kzz), '');
      h+=msbRow('k<sub>yz</sub> = k<sub>zz</sub>', '', msbR(B.kyz), '');
    }
    h+=msbRow(B.useB1? ('k<sub>zy</sub> = '+(B.c12? '0.6' : '0.8')+'k<sub>yy</sub>') : 'k<sub>zy</sub> = 1 &minus; 0.1&lambda;&#772;<sub>z</sub>U<sub>N.z</sub>/(C<sub>mLT</sub>&minus;0.25)', B.useB1? (B.c12? '0.6' : '0.8')+' x '+msbR(B.kyy) : msbEsc(B.kzyLbl)+'; &lambda;&#772;<sub>z</sub> = '+msbR(B.lamZ)+', U<sub>N.z</sub> = '+msbR(B.nz)+', C<sub>mLT</sub> = '+msbR(B.CmLT), msbR(B.kzy), B.useB1? 'Table B.1' : 'Table B.2');
    h+=msbRow('U<sub>Ny</sub>+k<sub>yy</sub>.U<sub>M.y</sub>+k<sub>yz</sub>.U<sub>M.z</sub>', msbR(B.ny)+'+'+msbR(B.kyy)+'x'+msbR(UMy)+'+'+msbR(B.kyz)+'x'+msbR(B.mzTerm), msbR(B.u1), nTag(B.u1));
    h+=msbRow('U<sub>Nz</sub>+k<sub>zy</sub>.U<sub>M.y</sub>+k<sub>zz</sub>.U<sub>M.z</sub>', msbR(B.nz)+'+'+msbR(B.kzy)+'x'+msbR(UMy)+'+'+msbR(B.kzz)+'x'+msbR(B.mzTerm), msbR(B.u2), nTag(B.u2));
  }

  /* ---- 5.8 Torsion Design ---- */
  if(T){
    h+=msbHead('Torsion Design');
    const tp=sec.tp||null;
    const Jv= tp&&tp.IT!=null? tp.IT : sec.J;
    const Hv= tp&&tp.Iw!=null? tp.Iw : sec.Iw;
    const av= T.p385? T.aa : (tp&&tp.a!=null? tp.a*1000 : null);
    h+=msbRow('J, H, a, Q<sub>f</sub>, Q<sub>w</sub>', g(Jv,2)+' cm&#8308;, '+msbDash(Hv,v=>g(v,4))+' dm&#8310;, '+msbDash(av,v=>g(v,0))+' mm, &mdash;, &mdash;', '', 'P385 App A');
    h+=msbRow('W<sub>n0</sub>, S<sub>w1</sub>', (tp&&tp.Wn0!=null? g(tp.Wn0,1)+' cm&sup2;' : '&mdash;')+', '+(tp&&tp.Sw1!=null? g(tp.Sw1,0)+' cm&#8308;' : '&mdash;')+(tp&&tp.Wn2!=null? '; W<sub>n2</sub> = '+g(tp.Wn2,1)+' cm&sup2;, S<sub>w2</sub> = '+g(tp.Sw2,1)+', S<sub>w3</sub> = '+g(tp.Sw3,1)+' cm&#8308;, e<sub>0</sub> = '+g(tp.e0,1)+', e<sub>sc</sub> = '+g(tp.esc,1)+' mm' : ''), '', 'P385 App A');
    const xT= T.p385&&T.cross? T.cross.x/1000 : (a.tors? a.tors.Tpos : 0);
    h+=msbSub('Torsion Bending Design @ '+msbM(xT)+' m');
    h+=msbRow('T<sub>Ed</sub> (max)', 'combination '+msbEsc(T.governT), msbKNm(T.TEd)+' kN.m', '');
    if(T.p385){
      // [beam-v03 addition, 19 Sep 2026 G4] method line: P385 closed forms where
      // they apply, the warping-torsion FE elsewhere (with its mesh error)
      if(T.fe){
        h+=msbRow('Torsion analysis', 'EI<sub>w</sub>&phi;&#8279; &minus; GI<sub>T</sub>&phi;&Prime; = m<sub>t</sub>(x): '+msbEsc(T.methodLabel)+'; '+T.bcText+'; closed forms not applicable: '+msbEsc((T.feReasons||[]).join('; '))+(T.nSolves!=null? '; '+T.nSolves+' FE solve(s)'+(T.nCached? ' + '+T.nCached+' cached' : '') : ''),
          'mesh error '+(T.meshError*100).toFixed(3)+' %', T.meshConverged? 'FE (&le; '+(T.meshBlock*100).toFixed(1)+' %)' : '<span class="ms-warn">BLOCKED</span>');
      } else {
        h+=msbRow('Torsion analysis', msbEsc(T.methodLabel)+'; '+T.bcText, 'L/a = '+f1(T.X,2), 'P385 App C');
      }
      h+=msbRow('&phi;<sub>max</sub> (ULS)', T.fe? 'warping-torsion FE, '+T.bcText : 'L/a = '+f1(T.X,2)+'; fork ends, warping free (P385 Cases 3/4/10)', f1(T.phiUmax,4)+' rad = '+f1(T.phiUmax*180/Math.PI,2)+'&deg;', '');
      h+=msbRow('B<sub>Ed</sub> = EI<sub>w</sub>&phi;&Prime; (max)', '@ x = '+msbM((T.BMaxPos||0)/1000)+' m', f1(T.BMax,3)+' kN.m&sup2;', '');
      h+=msbRow('M<sub>w.Ed</sub> = EI<sub>w</sub>&phi;&Prime;/(h&minus;t<sub>f</sub>)', 'max over span', msbKNm(T.MwMax)+' kN.m', '');
      // 20 Sep 2026 torsion + N/Mz: with an imposed M_z the line becomes M_z,tot = M_z,Ed + phi.M_y (imposed constant
      // moment + the largest coincident second-order minor-axis moment of the twist, SCI P385 3.1.2 / EN 1993-1-1
      // 5.2.1(3)); the station's M_y and phi are looked up in c.tor.grids (msbMzTwistRow), nothing recomputed
      const mzImp=(T.MzImp||0)>1e-9;
      if(mzImp){
        const rz=msbMzTwistRow(T);
        h+=msbRow('M<sub>z,tot</sub> = M<sub>z.Ed</sub> + &phi;.M<sub>y.Ed</sub>', msbKNm(T.MzImp)+' + '+(rz? msbKNm(rz.My)+' x '+g(Math.abs(rz.phi)*1000,2)+'e-3 (@ x = '+msbM(rz.x/1000)+' m)' : msbKNm(T.MzMax))+' = '+msbKNm(T.MzImp)+' + '+msbKNm(T.MzMax)+' (imposed constant M<sub>z.Ed</sub> + max coincident twist moment)', msbKNm(T.MzTot)+' kN.m', 'P385 3.1.2');
      } else h+=msbRow('M<sub>z.Ed</sub> = &phi;.M<sub>y.Ed</sub>', 'max coincident', msbKNm(T.MzMax)+' kN.m', '');
      const cr=T.cross;
      // 20 Sep 2026 torsion + N/Mz: the 3.1.2 M_z term is M_z,tot (cross.MzTot; = phi.M_y when no M_z is imposed)
      const crMz=(cr.MzTot!=null? cr.MzTot : cr.Mz), mzLbl=mzImp? 'M<sub>z,tot</sub>' : 'M<sub>z</sub>';
      h+=msbRow(T.cls12? '(M<sub>y</sub>/M<sub>pl.y</sub>)&sup2; + M<sub>w</sub>/M<sub>pl.f</sub> + '+mzLbl+'/M<sub>pl.z</sub>' : 'M<sub>y</sub>/M<sub>el.y</sub> + '+mzLbl+'/M<sub>el.z</sub> + M<sub>w</sub>/M<sub>f.Rd</sub>',
        '@ x = '+msbM(cr.x/1000)+' m: '+msbKNm(cr.My)+', '+msbKNm(cr.Mw)+', '+msbKNm(crMz)+' kN.m'+(mzImp? ' ('+mzLbl+' = '+msbKNm(cr.MzImp)+' + '+msbKNm(cr.Mz)+')' : '')+'; '+(T.cls12? 'M<sub>pl.y</sub> = '+msbKNm(T.Mply)+', M<sub>pl.f</sub> = '+msbKNm(T.Mplf)+', M<sub>pl.z</sub> = '+msbKNm(T.Mplz) : 'M<sub>el.y</sub> = '+msbKNm(T.Mely)+', M<sub>el.z</sub> = '+msbKNm(T.Melz)+', M<sub>f.Rd</sub> = '+msbKNm(T.Melf)),
        msbR(cr.u), msbWarn(cr.u<=1.0001)+' P385 3.1.2');
      // reference software "Combined Torsion buckling": the EN 1993-6 Annex A interaction with its amplifier k = kw.kzw.k_alpha (c.annex)
      if(AN){
        const unb = !isFinite(AN.kAlpha);
        // 20 Sep 2026 torsion + N/Mz: M_z,Ed of (A.1) = M_z,tot = imposed M_z + phi.M_y at the station, in the C_mz term and in
        // k_zw = 1 - M_z,tot/M_z,Rk (annexAEval); C_mz = 1.0 with an imposed M_z (constant diagram, psi = 1) - its reason printed
        const anImp=(AN.MzImp||0)>1e-9;
        if(anImp) h+=msbRow('C<sub>mz</sub> (EN 1993-6 A.1)', msbEsc(AN.CmzBasis||''), msbR(AN.Cmz), 'Table B.3');
        h+=msbRow('k = k<sub>w</sub>.k<sub>zw</sub>.k<sub>&alpha;</sub>', msbR(AN.kw)+' x '+msbR(AN.kzw)+' x '+(unb? '&infin; (M<sub>y.Ed</sub> &ge; M<sub>cr</sub> = '+msbKNm(AN.McrA)+' kN.m)' : msbR(AN.kAlpha))+(anImp&&!unb? ' (k<sub>zw</sub> = 1 &minus; M<sub>z,tot</sub>/M<sub>z.Rk</sub> = 1 &minus; '+msbKNm(AN.MzTot)+'/'+msbKNm(AN.MzR)+')' : ''), unb? '&mdash;' : msbR(AN.kw*AN.kzw*AN.kAlpha), 'EN 1993-6 A');
        h+=msbRow('M<sub>y</sub>/M<sub>b.Rd</sub> + C<sub>mz</sub>'+(anImp? 'M<sub>z,tot</sub>' : 'M<sub>z</sub>')+'/M<sub>z.Rk</sub> + k<sub>w</sub>k<sub>zw</sub>k<sub>&alpha;</sub>M<sub>w</sub>/M<sub>f.Rk</sub>',
          unb? 'k<sub>&alpha;</sub> unbounded: M<sub>y.Ed</sub> reaches M<sub>cr</sub> = '+msbKNm(AN.McrA)+' kN.m' : '@ x = '+msbM((AN.x||0)/1000)+' m: '+msbKNm(AN.My)+'/'+msbKNm(AN.MbA)+' + '+msbR(AN.Cmz)+'x'+msbKNm(AN.Mz)+'/'+msbKNm(AN.MzR)+' + '+msbR(AN.kw)+'x'+msbR(AN.kzw)+'x'+msbR(AN.kAlpha)+'x'+msbKNm(AN.Mw)+'/'+msbKNm(AN.MfR)+(anImp? ' (M<sub>z,tot</sub> = '+msbKNm(AN.MzImp)+' imposed + '+msbKNm(AN.MzTwist)+' twist &phi;.M<sub>y</sub>)' : ''),
          unb? '&ge; 99' : msbR(AN.u), unb? '<span class="ms-warn">Warning</span> FAIL: M<sub>y.Ed</sub> &ge; M<sub>cr</sub>' : msbWarn(AN.u<=1.0001)+' EN 1993-6 A');   // 20 Sep 2026: k_alpha unbounded is a FAIL (utilisation 99), not a blocked check
        // 20 Sep 2026 review: EN 1993-6 Annex A states (A.1) for I-section members; the channel application is the SCI P385 basis
        if(T.chan) h+=msbRow('(A.1) scope', 'EN 1993-6 Annex A (informative) states (A.1) for members with an I section; applied to the channel on the SCI P385 basis (P385 6.2/8.2 form with M<sub>f.Rd</sub> = t<sub>f</sub>b&sup2;f<sub>y</sub>/4) [verify]', '', 'EN 1993-6 A', 'ms-basis');
      }
      // 20 Sep 2026 torsion + N/Mz: the basis of the combined verification (tor.combinedBasis, engine text) with N_Ed
      // and / or an imposed M_z, then - N_Ed > 0 - the superposition of Eq 6.62 and (A.1) as the engine labels it:
      // information only, no OK / Warning tag, never a utilisation (tor.superposition)
      if(T.combinedActive && T.combinedBasis) h+=msbRow('Basis with N<sub>Ed</sub> / M<sub>z.Ed</sub>', msbEsc(T.combinedBasis), '', 'EN 1993-1-1, EN 1993-6', 'ms-basis');
      if(T.superposition){
        const sp=T.superposition, at=sp.at;
        h+=msbRow(msbEsc(sp.label), 'N<sub>Ed</sub>/(&chi;<sub>z</sub>N<sub>Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>zy</sub>M<sub>y.Ed</sub>/(&chi;<sub>LT</sub>M<sub>y.Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>zz</sub>M<sub>z,tot</sub>/(M<sub>z.Rk</sub>/&gamma;<sub>M1</sub>) + k<sub>w</sub>k<sub>zw</sub>k<sub>&alpha;</sub>M<sub>w</sub>/M<sub>f.Rd</sub> = '+msbR(sp.u62)+' (Eq 6.62'+(B? ' with M<sub>z,tot</sub> = '+msbKNm(B.MzEd) : '')+') + '+msbR(sp.uw)+' (warping term of (A.1)'+(at? ' @ x = '+msbM(at.x/1000)+' m: k<sub>w</sub> = '+msbR(at.kw)+', k<sub>zw</sub> = '+msbR(at.kzw)+', k<sub>&alpha;</sub> = '+msbR(sp.kAlpha)+(AN? '' : ' (no LTB check: M<sub>cr</sub> unbounded)')+', M<sub>w</sub> = '+msbKNm(at.Mw)+' kN.m' : '')+') =', msbR(sp.u), 'information only', 'ms-advrow');
      }
      h+=msbRow('End torques T<sub>t</sub>', 'St Venant part GI<sub>T</sub>&phi;&prime; at x = 0 / x = L'+(T.TEnds? '; total T = GI<sub>T</sub>&phi;&prime; &minus; EI<sub>w</sub>&phi;&#8244; = '+msbKNm(Math.abs(T.TEnds[0]))+' / '+msbKNm(Math.abs(T.TEnds[1]))+' kN.m' : ''), msbKNm(Math.abs(T.TtEnds[0]))+' / '+msbKNm(Math.abs(T.TtEnds[1]))+' kN.m', '');
      // Torsion Shear Design @ x (reference software sub-block): St Venant (+ warping, channel) shear stress at the governing V-T station, the shear-torsion reduction and the ratio
      const vt=T.vt||{};
      h+=msbSub('Torsion Shear Design @ '+msbM((vt.x||0)/1000)+' m');
      h+=msbRow('&tau;<sub>t</sub>'+(T.chan? ', &tau;<sub>w</sub>' : '')+' at the V-T station', '@ x = '+msbM((vt.x||0)/1000)+' m'+(vt.combo? ' ('+msbEsc(vt.combo)+')' : '')+': T<sub>t</sub> = '+msbKNm(vt.T!=null? vt.T : 0)+' kN.m, V = '+msbKN(vt.V!=null? vt.V : c.Fv)+' kN', f1(vt.tauT||0,2)+(T.chan? ' / '+f1(vt.tauW||0,2) : '')+' N/mm&sup2;', 'P385');
      h+=msbRow(T.chan? 'V<sub>pl.T.Rd</sub> = [&radic;(1 &minus; &tau;<sub>t</sub>/(1.25f<sub>y</sub>/&radic;3)) &minus; &tau;<sub>w</sub>/(f<sub>y</sub>/&radic;3)].V<sub>pl.Rd</sub>' : 'V<sub>pl.T.Rd</sub> = &radic;(1 &minus; &tau;<sub>t</sub>/(1.25f<sub>y</sub>/&radic;3)).V<sub>pl.Rd</sub>',
        '@ x = '+msbM((vt.x||0)/1000)+' m: &tau;<sub>t</sub> = '+f1(vt.tauT||0,2)+(T.chan? ' (&tau;<sub>w</sub> = '+f1(vt.tauW||0,2)+')' : '')+' N/mm&sup2;; V<sub>pl.Rd</sub> = '+msbKN(c.VcRd)+'; S<sub>mod</sub> = V<sub>pl.T.Rd</sub>/V<sub>pl.Rd</sub> = '+msbR(T.VplTRd/Math.max(c.VcRd,1e-9)), msbKN(T.VplTRd)+' kN', '6.2.7(9)');
      h+=msbRow('V<sub>Ed</sub>/V<sub>pl.T.Rd</sub>', msbKN(vt.V!=null? vt.V : c.Fv)+' / '+msbKN(T.VplTRd), T.vtZeroCapacity? '&infin;' : msbR(T.vtUtil), msbWarn(!T.vtZeroCapacity && T.vtUtil<=1.0001));
    } else if(T.box){
      // reference software box form: J, C (= W_t), tau_t.Ed = T/C, the torsion-modified local capacity (not evaluated here), then the torsion shear in the web
      h+=msbRow('W<sub>t</sub> (= C)', msbEsc(T.WtSrc)+'; I<sub>t</sub> = '+g(T.ItShow/1e4,1)+' cm&#8308;', g(T.Wt/1e3,1)+' cm&sup3;', '');
      h+=msbRow('T<sub>Rd</sub> = f<sub>y</sub>W<sub>t</sub>/(&radic;3&gamma;<sub>M0</sub>)', msbInt(c.fy)+' x '+g(T.Wt/1e3,1)+'/(&radic;3 x 1)', msbKNm(T.TRd)+' kN.m', '6.2.7(7)');
      h+=msbRow('T<sub>Ed</sub>/T<sub>Rd</sub>', msbKNm(T.TEd)+' / '+msbKNm(T.TRd), msbR(T.torUtil), msbWarn(T.torUtil<=1.0001));
      if(T.tauMax!=null) h+=msbRow('&tau;<sub>t.Ed</sub> = T<sub>Ed</sub>/W<sub>t</sub>', msbKNm(T.TEd)+' x 10&sup3;/'+g(T.Wt/1e3,1), f1(T.tauMax,2)+' N/mm&sup2;', '6.2.7(7)');
      // 20 Sep 2026 review: the verdict basis of a hollow section is the code's plastic route (T_Ed/T_Rd, V_pl.T.Rd with 6.2.8(4) rho and
      // 6.2.9.1/6.2.10); the elastic (6.1) check binds only for a Class 3 section (elasticBindingPolicy) and is information otherwise
      h+=msbRow('Modified Local Capacity (M<sub>y.Ed</sub>/(M<sub>pl.y.Rd</sub>.S<sub>mod</sub>))<sup>&alpha;</sup> + (M<sub>z.Ed</sub>/(M<sub>pl.z.Rd</sub>.S<sub>mod</sub>))<sup>&beta;</sup>', 'n/a - not evaluated by beam-v03 (the cl 6.2.7 T<sub>Ed</sub>/T<sub>Rd</sub>, V<sub>pl.T.Rd</sub> with the 6.2.8(4) &rho; in 6.2.9.1/6.2.10'+(T.elastic&&T.elastic.binding? ' and the elastic (6.1) check (Class 3) are' : ' checks are')+' the verdict basis'+(T.elastic&&!T.elastic.binding? '; the elastic (6.1) check is printed for information' : '')+')', '&mdash;', 'not evaluated');
      // 20 Sep 2026 torsion + N/Mz: with N_Ed and / or an imposed M_z the hollow section prints M_z,tot (St Venant twist) and the basis text
      // (20 Sep 2026 review: phi.M_y is the SCI P385 second-order term admitted through 5.2.1(3), not a 6.2.7(5) quantity - tagged so)
      if(T.combinedActive){
        if((T.MzImp||0)>1e-9) h+=msbRow('M<sub>z,tot</sub> = M<sub>z.Ed</sub> + &phi;.M<sub>y.Ed</sub>', msbKNm(T.MzImp)+' + '+msbKNm(T.MzTwistMax||0)+' (imposed constant M<sub>z.Ed</sub> + max coincident twist moment of the St Venant solution)', msbKNm(T.MzTot)+' kN.m', 'P385 3.1.2 / 5.2.1(3)');
        if(T.combinedBasis) h+=msbRow('Basis with N<sub>Ed</sub> / M<sub>z.Ed</sub>', msbEsc(T.combinedBasis), '', 'EN 1993-1-1, EN 1993-6', 'ms-basis');
      }
      const vt=T.vt||{};
      h+=msbSub('Torsion Shear Design @ '+msbM((vt.x||0)/1000)+' m');
      h+=msbRow('V<sub>pl.T.Rd</sub> = [1 &minus; &tau;<sub>t</sub>/(f<sub>y</sub>/&radic;3)].V<sub>pl.Rd</sub>', '@ x = '+msbM((vt.x||0)/1000)+' m'+(vt.combo? ' ('+msbEsc(vt.combo)+')' : '')+': T = '+msbKNm(vt.T!=null? vt.T : 0)+' kN.m, &tau;<sub>t</sub> = '+f1(vt.tau||0,2)+' N/mm&sup2;; V<sub>pl.Rd</sub> = '+msbKN(c.VcRd)+'; S<sub>mod</sub> = V<sub>pl.T.Rd</sub>/V<sub>pl.Rd</sub> = '+msbR(T.VplTRd/Math.max(c.VcRd,1e-9)), msbKN(T.VplTRd)+' kN', '6.2.7(9) Eq 6.28');
      h+=msbRow('V<sub>Ed</sub>/V<sub>pl.T.Rd</sub>', msbKN(vt.V!=null? vt.V : c.Fv)+' / '+msbKN(T.VplTRd), T.vtZeroCapacity? '&infin;' : msbR(T.vtUtil), msbWarn(!T.vtZeroCapacity && T.vtUtil<=1.0001));
    } else {
      h+=msbRow('P385 warping analysis', 'not covered for this arrangement (see NOT VERIFIED)', '&mdash;', '<span class="ms-warn">BLOCKED</span>');
    }
    h+=msbNotVerifiedRows(nv.torsion);
  }

  /* ---- blocking messages that no block above claims ---- */
  h+=msbNotVerifiedRows(nv.general);

  /* ---- 5.9 Deflection Check ---- */
  h+=msbHead('Deflection Check - Load Case '+caseD);
  const dCant = a.deflection? !!a.deflection.cant : isCant;
  const dAbs = a.deflection && a.deflection.abs!=null ? a.deflection.abs : null;
  const dAbsGov = !!(a.deflection && a.deflection.absGoverns);
  // limit label: Span/divisor (L/divisorCant when an end is vertically free), with the absolute cap when entered
  const limLbl=(cant,div)=>(cant? 'Tip &delta; &le; L/' : 'In-span &delta; &le; Span/')+msbInt(div)+(dAbs!=null? ' (&le; '+msbMM(dAbs)+' mm)' : '');
  const limVals=(sg)=>msbMM(Math.abs(sg.dmax))+' &le; '+(sg.absGoverns? msbMM(sg.abs)+' mm (absolute limit governs; '+g(sg.span,0)+' / '+msbInt(sg.divisor)+' = '+msbMM(sg.limSpan)+' mm)' : g(sg.span,0)+' / '+msbInt(sg.divisor)+' = '+msbMM(sg.limit)+' mm'+(sg.abs!=null? ' (absolute limit '+msbMM(sg.abs)+' mm not governing)' : ''))+' @ x = '+msbM(sg.dpos/1000)+' m';
  const gseg = a.deflection || {dmax:c.dmax,dpos:dfl.dpos,span:c.span,divisor:c.divisor,limit:c.dlimit,limSpan:c.dlimit,abs:null,absGoverns:false};
  h+=msbRow(limLbl(dCant,c.divisor), limVals(gseg)+(dCant? ' (vertically free end: tip deflection relative to the held end; L/'+msbInt(c.divisor)+' per UK NA to EN 1993-1-1 Table NA.2 [verify], cantilever row)' : ''), msbMM(c.dmax)+' mm', msbWarn(c.defOk)+(dAbsGov? ' abs' : ''));
  // 20 Sep 2026: with torsion active reference software prints the SLS twist here - "Torq in Case n @ x: theta_max = .. rad = .. deg <= 2.00 deg";
  // the 2 degree limit is P385 guidance and stays an advisory in beam-v03 (it does not enter c.utils)
  if(T && (T.p385 || T.box)){
    const twRad= T.p385? T.phiSer : T.phiMax, twDeg= T.p385? T.phiSerDeg : T.phiDeg, twX= T.p385? T.phiSerPos : T.phiPos;
    const twN=msbCaseIndex({label:T.governTw},true,a);
    h+=msbRow('Torq in Case '+(twN!=null? twN : '?')+' @ '+msbM(twX)+' m: &theta;<sub>max</sub> &le; 2.00&deg;', f1(twRad,4)+' rad = '+f1(twDeg,2)+'&deg; &le; 2.00&deg; ('+msbEsc(T.governTw)+'; P385 guidance, advisory: does not enter the verdict)', f1(twDeg,2)+'&deg;', twDeg<=2? 'OK' : '<span class="ms-warn">&gt; 2&deg;</span> advisory');
  }

  /* ---- 5.10 Unity bar ---- */
  const findU=(re)=>{ const u=utils.find(u=>re.test(u.name)); return u? u.val : null; };
  const cells=[];
  const push=(name,val)=>{ if(val!=null && isFinite(val)) cells.push({name,val}); else if(val===undefined) return; else cells.push({name,val:null}); };
  const deflU=c.dlimit>0? c.dmax/c.dlimit : null;
  if(AX){
    // 20 Sep 2026 review: the cells carry what the block rows print - N_Ed/N_c.Rd when the Class-4 A_eff applies,
    // U_M.y against the same M_b.Rd (W_el.y/W_pl.y applied) as the printed U_M.y row, and an em dash (not 0.000)
    // for the cl 6.3.3 cells of a tension brief, where the interaction is not evaluated
    push(AX.aeff&&AX.aeff.active? 'N_Ed/N_(c.Rd)' : 'N_Ed/N_(pl.Rd)', AX.nUtil);
    push('Local', AX.mUtil);
    push('UNyz', B? Math.max(B.ny,B.nz) : null);
    push('UMyz', B? Math.max(B.Mx/Math.max(B.MbRdEff!=null? B.MbRdEff : B.MbRdI,1e-9), B.mzTerm) : (ltbChecked? c.ltbUtil : c.momUtil));
    push('Ax+M_6.61', B? B.u1 : null);
    push('Ax+M_6.62', B? B.u2 : null);
    push('Deflection', deflU);
    push('V/Vpl', c.shearUtil);
    push('MA/Mc', c.momUtil);
    if(!fullRest) push('LTB', c.ltbUtil);
  } else {
    push('MA/Mc', c.momUtil);
    push('M_(y.Ed)/M_(b.Rd)', fullRest? c.momUtil : c.ltbUtil);
    push('Deflection', deflU);
    push('V/Vpl', c.shearUtil);
  }
  const coexU=findU(/6\.2\.8|Pure shear failure/), torU=findU(/^Torsion|Bending\+torsion cross-section/), vtU=findU(/Shear\+torsion/), anU=findU(/LTB\+torsion/);
  const webU=findU(/^Web transverse force  F_Ed/), web72U=findU(/^Web transverse force \+ bending/);
  const mvnU=findU(/6\.2\.10/), tfU=findU(/6\.3\.1\.4/);
  // 20 Sep 2026 torsion + N/Mz: its own cell, every torsion case; 20 Sep 2026 review: from c.utils when verdict-binding, else from
  // c.info (labelled "(info)", outside Max - 6.2.7(5) is permissive and the Class 1/2 plastic route governs)
  const elU=findU(/^Elastic yield criterion \(6\.1\) with torsion/);
  const elI=(c.info||[]).find(u=>/^Elastic yield criterion \(6\.1\) with torsion/.test(u.name));
  if(coexU!=null) push('M-V', coexU);
  if(mvnU!=null) push('M-V-N', mvnU);
  if(tfU!=null) push('N_b.T', tfU);
  if(webU!=null) push('F/F_Rd', webU);
  if(web72U!=null) push('Web 7.2', web72U);
  if(torU!=null) push('Torsion', torU);
  if(elU!=null) push('Yield 6.1', elU); else if(elI && isFinite(elI.val)) cells.push({name:'Yield 6.1 (info)', val:elI.val, info:true});
  if(vtU!=null) push('V+T', vtU);
  if(anU!=null) push('LTB+T', anU);
  const maxU=msbMaxExclDeflection(utils);
  cells.push({name:'Max', val:maxU, max:true});
  h+='<div class="ms-unity"><div class="ms-unity-head">'+cells.map(x=>'<div'+(x.max? ' class="ms-max"':'')+'>'+x.name+'</div>').join('')+'</div>'+
     '<div class="ms-unity-vals">'+cells.map(x=>'<div class="'+(x.max? 'ms-max ':'')+(x.val!=null&&x.val>1.0001&&!x.info? 'ms-warn':'')+'">'+(x.val==null? '&mdash;' : msbR(x.val))+'</div>').join('')+'</div></div>';

  /* ---- verdict footer ---- */
  h+='<div class="ms-verdict '+(c.pass? 'ms-pass' : 'ms-failv')+'">'+verdict+(c.gov? ' &mdash; governing '+msbEsc(c.gov.name)+' = '+msbR(c.gov.val) : '')+'</div>';
  if(unsupported.length) h+='<div class="ms-footer">'+unsupported.map(m=>'<div><span class="ms-warn">NOT COVERED:</span> '+msbEsc(m)+'</div>').join('')+'</div>';
  // 20 Sep 2026 torsion + N/Mz: an engine note that already opens with "ADVISORY - " (the superposition text) is not prefixed twice
  // 20 Sep 2026: a "FAIL:" note (k_alpha unbounded) is printed as a failure line in the footer, the rest as advisories
  if(c.advisory && c.advisory.length) h+='<div class="ms-footer ms-adv">'+c.advisory.map(m=> /^FAIL:/.test(String(m))
    ? '<div><span class="ms-warn">FAIL:</span> '+msbEsc(String(m).replace(/^FAIL:\s*/,''))+'</div>'
    : '<div><b>ADVISORY:</b> '+msbEsc(m).replace(/^ADVISORY\s*-\s*/,'')+'</div>').join('')+'</div>';

  return '<div class="ms-brief'+(c.pass? '' : ' ms-fail')+'">'+h+'</div>';
}


/* ==== beam-v03 module 16 ==== */
/* =====================================================================
   mcr-eigen-patch.js
   Replaces the C1-lookup / SN003a / SN006a / PFC-kappa LTB machinery in
   BEAM DESIGN CODE with a direct finite-element eigenvalue solution for
   the elastic critical moment.

   HOW TO APPLY
   ------------
   Paste this entire file inside the app's existing <script> block,
   immediately BEFORE the final `wire(); ... recompute();` lines.
   It reassigns two function bindings and injects one UI panel. Nothing
   else in the file is touched, so the encoding of the existing source is
   preserved.

   Mcr METHOD SWITCH  (S.mcrMethod: 'eigen' | 'standard')
   ---------------------------------------------------------
   The original closed-form function checksEC3UnrestrainedSCI() from
   js/checks/eurocode-checks.js is captured BEFORE the reassignment as
   window.checksEC3UnrestrainedStandard and the patched function delegates
   to it when S.mcrMethod === 'standard'. That is the STANDARD method:
   C1 from the NCCI SN003a tables / SCI end-moment curve / Serna quarter-point
   expression, the SN003a closed-form Mcr with C2*zg where published, SN006a
   for cantilevers and the P385/P362 channel kappa chain - as reference software-
   type software does. The eigen method (default) remains the FE eigenvalue
   solution for the actual moment diagram, load heights, restraints and hinges.
   Both methods expose ltb.McrStandard (closed form for the same segment) and
   ltb.McrEigen (null in standard mode: the eigensolver is not run) so the
   report can print "Mcr eigen / Mcr standard". The helpers below are
   therefore the standard-method implementation and must STAY ALIVE:

       C1_END_MOMENT, SN006 / sn006C, sernaC1, c1FromPsi, computeC1, mcrEC3
                            (js/01-computation-engine.js)
       sn003aC1, c1Inputs, c1Segment, mcrClosedForm, mcrSN006aFor,
       mcrStandardFor, cmTableB3, annexB2, checksEC3UnrestrainedSCI,
       annexAEval, torsionSuperposition, torsionCombinedBasis,
       ELASTIC_TORSION_UTIL_NAME
                            (js/checks/eurocode-checks.js; the last four
                            since 20 Sep 2026 torsion + N/Mz)

   WHAT CHANGES IN THE NUMBERS (eigen method vs the standard method)
   -----------------------------------------------------------------
   * C1 is no longer an input. Mcr comes out of the eigenproblem directly
     and lamLT = sqrt(Wy*fy/Mcr).
   * The old `computeC1`/`sn003aC1` returned C1 = 1.0 for every combined
     end-moment + transverse-load case. The true value is 1.13 to 2.75.
     Expect Mb,Rd to RISE on those beams. This is a correction, not a
     relaxation, but re-run your worked examples before trusting it.
   * The simplified route (P362 Expn 6.55) is gone. It was the design
     basis; the Mcr route was only used to rescue it.
   * Load height enters through zg on every load, for every moment
     diagram - not just the two shapes SN003a tabulates C2 for.
   * The destabilising x1.2 switch no longer affects EC3 LTB (zg does).
     LE-factor no longer affects EC3 LTB (restraint positions do).
     Both still drive the strut check.
   * PFC: lamLT now comes from Mcr, not (L/iz)/kappa. See NOTE ON CHANNELS.
   * Cantilevers: no SN006a table-range blocking (kwt <= 1, -2 <= eta <= 3).

   SIGN CONVENTIONS (must hold, and are asserted at load time)
   -----------------------------------------------------------
     x  along the span, 0..L        [mm]
     z  vertical, POSITIVE UPWARD, from the shear centre
     M  sagging POSITIVE            [N.mm]
     q,P transverse loads POSITIVE DOWNWARD
     zg is read directly from the UI as height above the shear centre;
     POSITIVE when the load acts ABOVE the shear centre
     zj > 0 when the compression flange is the larger flange

   Note the app's own comboLoads() returns downward loads as NEGATIVE.
   This patch reads S.loads directly and re-signs them; it does not use
   comboLoads(). assertSaggingPositive() below verifies the BMD sign at
   load time and throws loudly if a future edit flips it.

   NOTE ON CHANNELS
   ----------------
   A PFC bent about its major axis is symmetric about that axis: the
   mirror across mid-height maps top flange onto bottom flange. Hence
   zj = 0 and the shear centre lies on the horizontal centroidal axis.
   Its offset e0 is HORIZONTAL, and an eccentric load applies a PRIMARY
   TORQUE - a separate action, not a buckling term. It must not be folded
   into zg. So: eigen Mcr with zj = 0 and Iw about the shear centre, then
   the EN 1993-6 Annex A interaction (already in this app) combines it
   with the P385 torsion. If a channel carries eccentric load and the
   Annex A check cannot run, PASS is blocked.
   ===================================================================== */
(function () {
  'use strict';

  var PI = Math.PI, G_STEEL = 81000;   // N/mm2, per SN003a / P385 (NOT E/2.6)
  // Tool thresholds, not code requirements. Near lamLT ~ 1, chi_LT varies
  // roughly with sqrt(Mcr), so d(MbRd)/MbRd is about 0.5*d(Mcr)/Mcr.
  var MESH_WARN = 0.001;   // 0.1% - mention it
  var MESH_BLOCK = 0.005;  // 0.5% - refuse to certify

  /* ================================================================
     PART 1 - LTB eigenvalue engine
     Verified against: uniform moment (exact), Iw=0 (exact),
     Kitipornchai & Trahair monosymmetric closed form (exact),
     zg reversal identity (exact), midspan restraint -> L/2 problem,
     Timoshenko cantilever constants 4.013 / 12.85.
     ================================================================ */

  function zeros(n, m) { var A = [], i; for (i = 0; i < n; i++) A.push(new Float64Array(m === undefined ? n : m)); return A; }
  function transpose(A) { var n = A.length, m = A[0].length, B = zeros(m, n), i, j; for (i = 0; i < n; i++) for (j = 0; j < m; j++) B[j][i] = A[i][j]; return B; }
  function nrm2(x) { var s = 0, i; for (i = 0; i < x.length; i++) s += x[i] * x[i]; return Math.sqrt(s); }
  function matvec(A, x, out) { var n = A.length, i, j, s, Ai; for (i = 0; i < n; i++) { s = 0; Ai = A[i]; for (j = 0; j < n; j++) s += Ai[j] * x[j]; out[i] = s; } return out; }
  function seedVec(n, s) { var x = new Float64Array(n), i; for (i = 0; i < n; i++) x[i] = Math.sin(1.7 * (i + 1) + s) + 0.31 * Math.cos(0.37 * (i + 1) + s); var d = nrm2(x); for (i = 0; i < n; i++) x[i] /= d; return x; }

  var GP = [[0.033765242898424, 0.085662246189585], [0.169395306766868, 0.180380786524069],
            [0.380690406958402, 0.233956967286346], [0.619309593041598, 0.233956967286346],
            [0.830604693233132, 0.180380786524069], [0.966234757101576, 0.085662246189585]];

  function shape(xi, Le) {
    var x2 = xi * xi, x3 = x2 * xi;
    return {
      N:   [1 - 3 * x2 + 2 * x3, Le * (xi - 2 * x2 + x3), 3 * x2 - 2 * x3, Le * (-x2 + x3)],
      Np:  [(-6 * xi + 6 * x2) / Le, 1 - 4 * xi + 3 * x2, (6 * xi - 6 * x2) / Le, -2 * xi + 3 * x2],
      Npp: [(-6 + 12 * xi) / (Le * Le), (-4 + 6 * xi) / Le, (6 - 12 * xi) / (Le * Le), (-2 + 6 * xi) / Le]
    };
  }
  function kBend(EI, Le) {
    var c = EI / (Le * Le * Le), L = Le, L2 = Le * Le;
    return [[12 * c, 6 * L * c, -12 * c, 6 * L * c], [6 * L * c, 4 * L2 * c, -6 * L * c, 2 * L2 * c],
            [-12 * c, -6 * L * c, 12 * c, -6 * L * c], [6 * L * c, 2 * L2 * c, -6 * L * c, 4 * L2 * c]];
  }
  function kTors(GJ, Le) {
    var c = GJ / (30 * Le), L = Le, L2 = Le * Le;
    return [[36 * c, 3 * L * c, -36 * c, 3 * L * c], [3 * L * c, 4 * L2 * c, -3 * L * c, -L2 * c],
            [-36 * c, -3 * L * c, 36 * c, -3 * L * c], [3 * L * c, -L2 * c, -3 * L * c, 4 * L2 * c]];
  }
  function cholesky(A) {
    var n = A.length, L = zeros(n), i, j, k, s;
    for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
      s = A[i][j]; for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (s <= 1e-300) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
    }
    return L;
  }
  function solveLower(L, B) {
    var n = L.length, m = B[0].length, X = zeros(n, m), i, j, c, s;
    for (c = 0; c < m; c++) for (i = 0; i < n; i++) {
      s = B[i][c]; for (j = 0; j < i; j++) s -= L[i][j] * X[j][c]; X[i][c] = s / L[i][i];
    }
    return X;
  }
  function solveUpperT(L, b) {
    var n = L.length, x = new Float64Array(n), i, j, s;
    for (i = n - 1; i >= 0; i--) { s = b[i]; for (j = i + 1; j < n; j++) s -= L[j][i] * x[j]; x[i] = s / L[i][i]; }
    return x;
  }

  /* Spectral radius via power iteration on A^2.
     Necessary because when zg = zj = 0 the elastic matrix has no v-phi
     coupling, Cholesky preserves the (v,phi) block split, Kg is purely
     off-diagonal, and A = [[0,C],[C^T,0]] has an exactly symmetric +-mu
     spectrum on which plain power iteration silently returns a wrong,
     "converged" Rayleigh quotient. Squaring makes the dominant eigenvalue
     simple; shifting afterwards keeps it that way. */
  function spectralRadius(A) {
    var n = A.length, x = seedVec(n, 0.3), v = new Float64Array(n), w = new Float64Array(n);
    var r2 = 0, r2o = NaN, i, it, nw;
    for (it = 0; it < 5000; it++) {
      matvec(A, x, v);
      r2 = 0; for (i = 0; i < n; i++) r2 += v[i] * v[i];
      matvec(A, v, w); nw = nrm2(w);
      if (nw < 1e-300 || r2 < 1e-300) return Math.sqrt(Math.max(r2, 0));
      for (i = 0; i < n; i++) x[i] = w[i] / nw;
      if (it > 3 && Math.abs(r2 - r2o) <= 1e-15 * r2) break;
      r2o = r2;
    }
    return Math.sqrt(r2);
  }
  function dominantShifted(A, shift, seed) {
    var n = A.length, x = seedVec(n, seed), y = new Float64Array(n), Ax = new Float64Array(n);
    var mu = 0, muOld = NaN, i, it, ny, conv = false;
    for (it = 0; it < 20000; it++) {
      matvec(A, x, Ax);
      for (i = 0; i < n; i++) y[i] = Ax[i] - shift * x[i];
      ny = nrm2(y); if (ny < 1e-300) { mu = shift; conv = true; break; }
      for (i = 0; i < n; i++) x[i] = y[i] / ny;
      matvec(A, x, Ax);
      mu = 0; for (i = 0; i < n; i++) mu += x[i] * Ax[i];
      if (it > 3 && Math.abs(mu - muOld) <= 1e-14) { conv = true; break; }
      muOld = mu;
    }
    return { mu: mu, vec: x, converged: conv };
  }

  function meshNodes(L, forced, nElem) {
    var pts = [0, L], i, j;
    forced.forEach(function (x) { if (x > 1e-9 && x < L - 1e-9) pts.push(x); });
    pts.sort(function (a, b) { return a - b; });
    var uniq = [pts[0]];
    for (i = 1; i < pts.length; i++) if (pts[i] - uniq[uniq.length - 1] > 1e-6) uniq.push(pts[i]);
    var nodes = [uniq[0]];
    for (i = 0; i < uniq.length - 1; i++) {
      var a = uniq[i], b = uniq[i + 1], k = Math.max(1, Math.round(nElem * (b - a) / L));
      for (j = 1; j <= k; j++) nodes.push(a + (b - a) * j / k);
    }
    return nodes;
  }
  function nearestNode(nodes, x) {
    var bi = 0, bd = Infinity, i;
    for (i = 0; i < nodes.length; i++) { var d = Math.abs(nodes[i] - x); if (d < bd) { bd = d; bi = i; } }
    return bi;
  }

  function mcrOnce(p) {
    var E = p.E, G = p.G === undefined ? G_STEEL : p.G;
    var Iz = p.Iz, It = p.It, Iw = p.Iw || 0, L = p.L, zj = p.zj || 0;
    var nElem = p.nElem || 32, dl = p.distLoads || [], pl = p.pointLoads || [], rst = p.restraints || [];
    var M = p.moment, i, j, k;

    var forced = [];
    rst.forEach(function (r) { forced.push(r.x); });
    pl.forEach(function (q) { forced.push(q.x); });
    dl.forEach(function (q) { forced.push(q.x1); forced.push(q.x2); });
    /* an applied couple makes M(x) jump; a node there keeps the jump on an
       element boundary, where the Hermite O(h^4) extrapolation holds (19 Sep
       2026 review finding F-D: a couple off a node inflated the mesh error) */
    (p.couplePositions || []).forEach(function (x) { forced.push(x); });
    var nodes = meshNodes(L, forced, nElem), nn = nodes.length, nd = 4 * nn;

    function qzg(x) {
      var s = 0, d, w;
      for (var t = 0; t < dl.length; t++) {
        d = dl[t];
        if (x < d.x1 - 1e-9 || x > d.x2 + 1e-9) continue;
        w = (d.x2 - d.x1 < 1e-9) ? d.w1 : d.w1 + (d.w2 - d.w1) * (x - d.x1) / (d.x2 - d.x1);
        s += w * (d.zg || 0);
      }
      return s;
    }

    var Ke = zeros(nd), Kg = zeros(nd), vL = [0, 1, 4, 5], pLoc = [2, 3, 6, 7];
    for (var e = 0; e < nn - 1; e++) {
      var x0 = nodes[e], Le = nodes[e + 1] - x0;
      var Kb = kBend(E * Iz, Le), Kw = kBend(E * Iw, Le), Kt = kTors(G * It, Le);
      var g = []; for (i = 0; i < 8; i++) g.push(4 * e + i);
      for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        Ke[g[vL[i]]][g[vL[j]]] += Kb[i][j];
        Ke[g[pLoc[i]]][g[pLoc[j]]] += Kw[i][j] + Kt[i][j];
      }
      var cpl = zeros(4, 4), gpp = zeros(4, 4);
      for (k = 0; k < GP.length; k++) {
        var xi = GP[k][0], wt = GP[k][1] * Le, x = x0 + xi * Le;
        var S2 = shape(xi, Le), Mx = M(x), qz = qzg(x);
        for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
          cpl[i][j] += Mx * S2.Npp[i] * S2.N[j] * wt;                       // INT M v'' phi
          gpp[i][j] += (2 * zj * Mx * S2.Np[i] * S2.Np[j]                   // Wagner
                        - qz * S2.N[i] * S2.N[j]) * wt;                     // load height
        }
      }
      for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        Kg[g[vL[i]]][g[pLoc[j]]] += cpl[i][j];
        Kg[g[pLoc[j]]][g[vL[i]]] += cpl[i][j];
        Kg[g[pLoc[i]]][g[pLoc[j]]] += gpp[i][j];
      }
    }
    pl.forEach(function (q) {
      var ni = nearestNode(nodes, q.x);
      Kg[4 * ni + 2][4 * ni + 2] += -q.P * (q.zg || 0);
    });

    var fixed = {};
    rst.forEach(function (r) {
      var ni = nearestNode(nodes, r.x);
      if (r.v === undefined ? true : !!r.v) fixed[4 * ni + 0] = 1;
      if (r.vp) fixed[4 * ni + 1] = 1;
      if (r.phi === undefined ? true : !!r.phi) fixed[4 * ni + 2] = 1;
      if (r.phip) fixed[4 * ni + 3] = 1;
    });
    var free = []; for (var d2 = 0; d2 < nd; d2++) if (!fixed[d2]) free.push(d2);
    var nf = free.length;
    if (nf < 2) throw new Error('LTB model has fewer than 2 free degrees of freedom.');

    var Kef = zeros(nf), Kgf = zeros(nf);
    for (i = 0; i < nf; i++) for (j = 0; j < nf; j++) { Kef[i][j] = Ke[free[i]][free[j]]; Kgf[i][j] = Kg[free[i]][free[j]]; }

    var Lc = cholesky(Kef);
    if (!Lc) throw new Error('LTB elastic stiffness is singular: the member is laterally under-restrained. ' +
      'Provide at least two points of lateral restraint, or lateral-bending restraint at a built-in support.');

    var X = solveLower(Lc, Kgf), Y = solveLower(Lc, transpose(X)), A = transpose(Y);
    for (i = 0; i < nf; i++) for (j = i + 1; j < nf; j++) { var m = 0.5 * (A[i][j] + A[j][i]); A[i][j] = m; A[j][i] = m; }

    var rho = spectralRadius(A);
    if (!(rho > 0)) throw new Error('LTB geometric stiffness is null: the moment diagram is identically zero.');
    var As = zeros(nf);
    for (i = 0; i < nf; i++) for (j = 0; j < nf; j++) As[i][j] = A[i][j] / rho;
    var sig = 1 + 1e-6;
    var dLo = dominantShifted(As, sig, 0.7), dHi = dominantShifted(As, -sig, 1.9);
    var muMin = dLo.mu * rho, muMax = dHi.mu * rho;

    var Mref = 0;
    for (i = 0; i <= 2000; i++) { var av = Math.abs(M(L * i / 2000)); if (av > Mref) Mref = av; }

    var eps = 1e-12 * Math.max(Math.abs(muMin), Math.abs(muMax), 1e-300);
    var lamPos = (muMin < -eps) ? -1 / muMin : Infinity;
    var lamNeg = (muMax > eps) ? -1 / muMax : -Infinity;

    var mode = null;
    if (dLo.vec) {
      var dv = solveUpperT(Lc, dLo.vec), full = new Float64Array(nd);
      for (i = 0; i < nf; i++) full[free[i]] = dv[i];
      var pk = 0; for (i = 0; i < nn; i++) pk = Math.max(pk, Math.abs(full[4 * i + 2]));
      mode = { x: nodes.slice(), phi: [], v: [] };
      for (i = 0; i < nn; i++) { mode.phi.push(pk > 0 ? full[4 * i + 2] / pk : 0); mode.v.push(full[4 * i]); }
    }
    return { Mcr: lamPos * Mref, lambda: lamPos, McrRev: Math.abs(lamNeg) * Mref, lambdaRev: lamNeg,
             Mref: Mref, nElem: nn - 1, converged: dLo.converged && dHi.converged, mode: mode };
  }

  /* Richardson extrapolation on the O(h^4) eigenvalue error of cubic Hermite. */
  function mcrEigen(p) {
    if (!p.refine) return mcrOnce(p);
    var n0 = p.nElem || 32;
    var co = mcrOnce(Object.assign({}, p, { nElem: n0, refine: false }));
    var fi = mcrOnce(Object.assign({}, p, { nElem: 2 * n0, refine: false }));
    function rich(a, b) { return (isFinite(a) && isFinite(b)) ? (16 * b - a) / 15 : b; }
    fi.lambda = rich(co.lambda, fi.lambda); fi.Mcr = fi.lambda * fi.Mref;
    fi.lambdaRev = rich(co.lambdaRev, fi.lambdaRev); fi.McrRev = Math.abs(fi.lambdaRev) * fi.Mref;
    fi.meshError = Math.abs(fi.Mcr - co.lambda * co.Mref) / Math.max(Math.abs(fi.Mcr), 1e-9);
    fi.converged = co.converged && fi.converged;
    fi.refined = true;
    return fi;
  }

  function momentFromSamples(xs, Ms) {
    return function (x) {
      var n = xs.length;
      if (x <= xs[0]) return Ms[0];
      if (x >= xs[n - 1]) return Ms[n - 1];
      var lo = 0, hi = n - 1;
      while (hi - lo > 1) { var mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
      var t = (x - xs[lo]) / (xs[hi] - xs[lo]);
      return Ms[lo] + (Ms[hi] - Ms[lo]) * t;
    };
  }

  /* ================================================================
     PART 2 - sign-convention assertion
     Runs once at load. If a future edit flips comboLoads() or sfdBmd(),
     this fails loudly instead of quietly inverting the zg term (which
     would make top-flange loading look STABILISING).
     ================================================================ */
  function assertSaggingPositive() {
    var L = 6000, EI = 210000 * 3.71e8;
    var sup = [{ pos: 0, type: 'pinned' }, { pos: L, type: 'pinned' }];
    var loads = [{ type: 'udl', x1: 0, x2: L, w1: -10, w2: -10 }];  // downward, app's sign
    var r = solveBeam(L, EI, sup, loads);
    var fb = sfdBmd(L, sup, loads, r.reactions);
    var Mmid = interpAt(fb.xs, fb.M, L / 2);
    if (!(Mmid > 0)) throw new Error(
      'mcr-eigen-patch: BMD sign convention check FAILED. A simply supported beam under downward ' +
      'UDL returned M(L/2) = ' + Mmid + ', expected sagging POSITIVE. The zg load-height term would ' +
      'be inverted. Fix the convention before using this patch.');
    return true;
  }

  /* ================================================================
     PART 3 - build the eigen model from the app's state
     ================================================================ */

  /* Section properties, mm units. Prefers the P385 / section tables values
     (tp.IT, tp.Iw) which are taken about the SHEAR CENTRE for channels. */
  function secProps(sec) {
    var tp = sec.tp || {};
    return {
      Iz: sec.Iy * 1e4,                                        // minor-axis I, mm4
      It: ((tp.IT != null ? tp.IT : sec.J) || 0) * 1e4,        // St Venant, mm4
      Iw: ((tp.Iw != null ? tp.Iw : sec.Iw) || 0) * 1e12,      // warping, mm6 (about shear centre)
      hs: sec.D - sec.tf                                       // flange centroid separation
    };
  }

  /* zj: zero for every section in this app's library.
     - I/H (UB, UC): doubly symmetric.
     - PFC: symmetric about the MAJOR axis (top flange maps onto bottom).
     - SHS/RHS: doubly symmetric.
     Non-zero zj belongs to unequal-flange plated I-sections and tees.
     S.zj is exposed so a future plated-section path can supply it. */
  function zjFor(sec) { return (S.zj != null && isFinite(S.zj)) ? +S.zj : 0; }

  function ltbRestraintsFor(a) {
    /* End boundary conditions straight from the end degree-of-freedom flags
       (19 Sep 2026 single-span scope): v from U_y, v' from R_z, phi from R_x,
       phi' from warping (ltbEndRestraints, js/03-state-ui.js). A free end
       (cantilever tip) contributes no restraint at x = L. Intermediate lateral
       restraints follow with their own four flags. A warping flag is applied
       only where the section has a warping constant (warpingApplies,
       js/03-state-ui.js): on an I_w = 0 box the twist equation is second
       order and phi' is not a boundary condition (finding F-E). */
    var warp = warpingApplies(a.sec);
    var out = ltbEndRestraints(S).map(function (r) { return { x: r.x, v: r.v, vp: r.vp, phi: r.phi, phip: warp ? r.phip : 0 }; });
    (S.ltbRestraints || []).forEach(function (r) {
      out.push({ x: (+r.pos) * 1000, v: r.v !== false, vp: !!r.vp, phi: r.phi !== false, phip: warp ? !!r.phip : false });
    });
    return out;
  }

  /* Transverse loads of ONE ULS combination, re-signed DOWNWARD
     POSITIVE, each carrying the load-height zg. Moment loads contribute to
     M(x) but have no load-height term, so they are excluded from the load
     lists; their positions are returned as `couples` so the eigen mesh puts a
     node at every jump of the moment diagram.
     This deliberately does NOT call comboLoads(), whose sign is inverted. */
  function unitLoadsFor(a, combo) {
    var fac = combo.factors;
    var dl = [], pl = [], couples = [];
    comboLoadPieces(combo).forEach(function (p) {
      var f = p.factor;
      if (!f) return;
      if (p.type === 'moment') { if (Math.abs(p.M) > 1e-12) couples.push(p.pos); return; }
      var zg = typeof loadZgValue === 'function' ? loadZgValue(p.ld) : (+S.za || 0);
      if (p.type === 'point') pl.push({ x: p.pos, P: p.P * f * 1000, zg: zg });
      else dl.push({ x1: p.x1, x2: p.x2, w1: p.w1 * f, w2: p.w2 * f, zg: zg });
    });
    var gF = fac.G != null ? fac.G : 0;
    var sw = selfWeightValue(a.sec);
    if (gF !== 0 && sw > 0) dl.push({ x1: 0, x2: a.L, w1: sw * gF, w2: sw * gF, zg: 0 }); // self-weight acts at the centroid
    return { distLoads: dl, pointLoads: pl, couples: couples };
  }

  /* Full LTB solve for ONE ULS combination (res = {combo, fb}). Returns Mcr,
     plus a GENERALISED C1 defined as
        C1 = Mcr(actual diagram, zg=0, zj=0) / Mcr(uniform moment, same restraints)
     which reduces to the textbook C1 on a fork-fork span and remains
     meaningful for intermediately restrained and non-fork-ended members,
     where no tabulated C1 exists. Used only for kc (NA 2.18).
     `shared` caches the combo-independent uniform-moment reference solve so a
     multi-combination run does not repeat it. */
  /* Memo of eigen solves (19 Sep 2026 review, performance): keyed by every
     input of the solve - section constants, L, restraints, mesh, zj, the load
     lists and the sampled moment diagram - so a re-render after an unrelated
     input (deflection divisor, grade, s_s) or an identical pattern of another
     parent combination does not repeat a coarse + fine eigen solve. Bounded
     (oldest entry dropped); `shared.stats` counts solves and cache hits. */
  var EIGEN_CACHE = new Map(), EIGEN_CACHE_MAX = 96;
  function mcrEigenMemo(p, momentKey, stats) {
    var key = JSON.stringify([p.E, p.G, p.Iz, p.It, p.Iw, p.L, p.restraints, p.nElem, !!p.refine, p.zj, p.distLoads || [], p.pointLoads || [], p.couplePositions || [], momentKey]);
    var hit = EIGEN_CACHE.get(key);
    if (hit) { EIGEN_CACHE.delete(key); EIGEN_CACHE.set(key, hit); if (stats) stats.cached++; return hit; }
    var sol = mcrEigen(p);
    if (stats) stats.solved++;
    EIGEN_CACHE.set(key, sol);
    if (EIGEN_CACHE.size > EIGEN_CACHE_MAX) EIGEN_CACHE.delete(EIGEN_CACHE.keys().next().value);
    return sol;
  }
  function solveLTB(a, sec, res, shared) {
    var sp = secProps(sec), zj = zjFor(sec), ul = unitLoadsFor(a, res.combo);
    var gfb = res.fb;
    var base = { E: a.E, G: G_STEEL, Iz: sp.Iz, It: sp.It, Iw: sp.Iw, L: a.L,
                 restraints: ltbRestraintsFor(a), nElem: 32, refine: true, couplePositions: ul.couples };
    var stats = shared && shared.stats;
    var mKey = gfb.xs.length + ':' + gfb.M.join(',');

    var actual = mcrEigenMemo(Object.assign({}, base, {
      moment: momentFromSamples(gfb.xs, gfb.M), zj: zj,
      distLoads: ul.distLoads, pointLoads: ul.pointLoads
    }), mKey, stats);

    var shapeOnly = mcrEigenMemo(Object.assign({}, base, {
      moment: momentFromSamples(gfb.xs, gfb.M), zj: 0,
      distLoads: ul.distLoads.map(function (d) { return Object.assign({}, d, { zg: 0 }); }),
      pointLoads: ul.pointLoads.map(function (q) { return Object.assign({}, q, { zg: 0 }); })
    }), mKey, stats);

    var uniform = (shared && shared.uniform) ||
      mcrEigenMemo(Object.assign({}, base, { moment: function () { return 1e6; }, zj: 0 }), 'uniform', stats);
    if (shared) shared.uniform = uniform;

    var C1 = uniform.Mcr > 0 ? shapeOnly.Mcr / uniform.Mcr : 1;
    var zgVals = ul.distLoads.map(function(d){ return d.zg || 0; }).concat(ul.pointLoads.map(function(q){ return q.zg || 0; }));
    var zgRep = 0;
    zgVals.forEach(function(z){ if(Math.abs(z)>Math.abs(zgRep)) zgRep = z; });
    var zgUnique = [];
    zgVals.forEach(function(z){
      if(!zgUnique.some(function(u){ return Math.abs(u-z)<1e-9; })) zgUnique.push(z);
    });
    return { Mcr: actual.Mcr, McrRev: actual.McrRev, mode: actual.mode,
             meshError: actual.meshError, mcrConverged: actual.converged,
             c1Converged: shapeOnly.converged && uniform.converged,
             C1: C1, McrUniform: uniform.Mcr, McrShape: shapeOnly.Mcr,
             zg: zgRep, zgValues: zgUnique, zgUniform: zgUnique.length <= 1, zj: zj, nElem: actual.nElem, sp: sp };
  }

  function ltbCurve(sec) {
    // UK NA Table NA.1 (cl 6.3.2.3): hot-finished hollow sections share the
    // I/H h/b allocation; cold-formed hollow sections use c (h/b<=2) or d;
    // channels d. ONE allocation for both Mcr methods: ltbCurveNA() in
    // js/checks/eurocode-checks.js (the standard route uses the same function).
    if (typeof ltbCurveNA === 'function') return ltbCurveNA(sec);
    if (sec.isBox && sec.boxType === 'CF') return sec.D/sec.B <= 2 ? { alphaLT: 0.49, curve: 'c' } : { alphaLT: 0.76, curve: 'd' };
    if (sec.kind === 'channel') return { alphaLT: 0.76, curve: 'd' }; // not doubly symmetric
    var hb = sec.D / sec.B;
    return hb <= 2 ? { alphaLT: 0.34, curve: 'b' } : hb <= 3.1 ? { alphaLT: 0.49, curve: 'c' } : { alphaLT: 0.76, curve: 'd' };
  }

  /* ================================================================
     PART 4 - replacement LTB check
     ================================================================ */
  /* The closed-form (standard) implementation from eurocode-checks.js is kept
     under its own name; the patched binding delegates to it on request. */
  var checksEC3UnrestrainedStandard = window.checksEC3UnrestrainedSCI;
  if (typeof checksEC3UnrestrainedStandard !== 'function')
    throw new Error('mcr-eigen-patch: checksEC3UnrestrainedSCI (the standard closed-form method) must be loaded before this patch.');
  window.checksEC3UnrestrainedStandard = checksEC3UnrestrainedStandard;
  function mcrMethod() { return (S.mcrMethod === 'standard') ? 'standard' : 'eigen'; }

  window.checksEC3UnrestrainedSCI = function (a) {
    if (mcrMethod() === 'standard') return checksEC3UnrestrainedStandard(a);
    var b = checksEC3Restrained(a);
    var sec = a.sec, fy = a.py, gM1 = 1.0, Wy = b.Wy;
    var unsupported = b.unsupported.slice();
    var isCant = isCantilever(S);
    var ltb, warn = [];

    /* Closed hollow sections are no longer blanket-exempted. EN 1993-1-1
       cl 6.3.2.1(2) genuinely exempts square/circular hollow sections, but a
       slender RECTANGULAR hollow section on a long unrestrained span has a
       finite Mcr. The eigen solver handles the closed section directly
       (Iw ~ 0, large It): a square/stocky box returns lamLT far below 0.4 and
       LTB is "ignored" (cl 6.3.2.2(4)), reproducing the exemption exactly
       where it is genuine; a slender RHS gets a real chi_LT. */

    /* LTB is checked for EVERY enabled ULS combination, not only the one with
       the largest |Mmax|: Mcr depends on the moment-diagram SHAPE, so a combo
       with a slightly smaller peak moment but a more uniform diagram (lower
       Mcr) can govern the LTB utilisation. The uniform-moment reference solve
       is combo-independent and shared. */
    var combosLTB = (a.ulsResults || []).filter(function (r) { return Math.abs(r.Mmax) > 1e-9; });
    if (!combosLTB.length) combosLTB = [a.governM];
    var shared = { uniform: null, stats: { solved: 0, cached: 0 } };
    var evals = [];
    try {
      combosLTB.forEach(function (res) {
        evals.push({ res: res, sol: solveLTB(a, sec, res, shared) });
      });
    } catch (err) {
      unsupported.push('Elastic critical moment: ' + err.message);
      var stdFail = null;
      try { stdFail = mcrStandardFor(a, sec); } catch (e3) { stdFail = null; }
      ltb = { eigen: true, failed: true, err: err.message, MbRd: 0, Mcr: 0, C1: 1, kc: 1,
              curve: ltbCurve(sec), ign: false, chi: 0, f: 1, chiMod: 0, lamLT: 0, warn: warn,
              mcrMethod: 'eigen', McrEigen: null, McrRatio: null,
              McrStandard: (stdFail && stdFail.Mcr != null) ? stdFail.Mcr : null,
              std: stdFail, c1in: stdFail ? stdFail.c1in : null, c1seg: stdFail ? stdFail.seg : null, c1route: stdFail ? stdFail.route : null };
      return Object.assign({}, b, { sci: false, sciU: true, mcrMethod: 'eigen', unsupported: unsupported, ltb: ltb,
        ltbUtil: 99, ltbBasis: 'Mcr could not be computed', C1: 1, c1label: 'n/a',
        LE: a.L, utils: [{ name: 'LTB', val: 99 }], gov: { name: 'LTB', val: 99 }, pass: false,
        annex: null, buck: null });
    }

    var meshErrorMax = 0, allConverged = true;
    evals.forEach(function (ev) {
      meshErrorMax = Math.max(meshErrorMax, ev.sol.meshError || 0);
      if (!ev.sol.mcrConverged) allConverged = false;
    });
    if (!allConverged) {
      unsupported.push('Elastic critical moment: the eigensolver did not reach convergence tolerance at one or both mesh levels' +
        (evals.length > 1 ? ' for at least one load combination' : '') + '. ' +
        'M<sub>cr</sub> is not reliable for design and PASS is blocked. Verify M<sub>cr</sub> independently.');
    }
    if (meshErrorMax > MESH_BLOCK) {
      unsupported.push('Elastic critical moment: mesh convergence error is ' + (meshErrorMax * 100).toFixed(2) +
        ' %, above the ' + (MESH_BLOCK * 100).toFixed(1) + ' % limit this tool will certify. PASS is blocked.');
    } else if (meshErrorMax > MESH_WARN) {
      warn.push('Mesh convergence error ' + (meshErrorMax * 100).toFixed(2) + ' % (below the blocking limit, but worth noting).');
    }

    var zgAny = evals.some(function (ev) { return Math.abs(ev.sol.zg || 0) > 1e-9; });
    if (zgAny && S.destab) warn.push('The destabilising x1.2 switch is ignored on the EC3 path: load height is carried exactly by per-load zg. Untick it to avoid confusion.');
    /* 19 Sep 2026 verification campaign (UB-40): the switch declares a destabilising load, but with every
       z_g = 0 the eigenvalue would take the load at the shear centre and print a PASS the closed-form route
       (L_E x 1.2) refuses - a contradictory input, so PASS is blocked until z_g is entered or the switch cleared. */
    if (!zgAny && S.destab) unsupported.push('Destabilising loading is ticked but every load height z<sub>g</sub> is 0 (load at the shear centre): the eigenvalue M<sub>cr</sub> carries the load height exactly through z<sub>g</sub> and does not apply the &times;1.2 L<sub>E</sub> device of the closed-form route, so this solve would treat the load as non-destabilising. Enter the load height (e.g. z<sub>g</sub> = +h/2 for a top-flange load, with eccentricity/height inputs on) or untick the switch; PASS is blocked.');
    if (S.leFactor != null && S.leFactor !== '' && Math.abs(+S.leFactor - 1) > 1e-9) warn.push('The L<sub>E</sub> factor does not affect EC3 LTB on the eigen route; the buckling length is set by the end restraints and the restraint positions. It overrides the strut lengths of both axes (the end-fixity defaults are printed with them).');
    var warpFlagged = endsList(S).filter(function (e) { return e.warp; }).map(function (e) { return 'End ' + e.n; })
      .concat((S.ltbRestraints || []).filter(function (r) { return r.phip; }).map(function (r) { return 'restraint at x = ' + g(+r.pos, 2) + ' m'; }));
    if (warpFlagged.length && !warpingApplies(sec)) warn.push('Warping flag at ' + warpFlagged.join(', ') + ' not applied: this closed section has I<sub>w</sub> = 0, so its twist equation is the second-order St Venant form and &phi;&prime; is not a boundary condition (EN 1993-1-1 6.2.7(7): warping of closed hollow sections may be neglected). The eigen model holds v, v&prime; and &phi; as flagged.');

    var curve = ltbCurve(sec);

    /* kc's C1 must follow the NA 2.18 basis: the moment shape BETWEEN
       RESTRAINTS. With intermediate restraints the whole-member shape ratio
       can exceed the critical bay's own C1 (e.g. a near-uniform critical bay
       inside a varied overall diagram), which would overstate the f-factor
       benefit. So each combination's kc uses the C1 of the bay its OWN
       eigenmode localises in - isolated fork-fork, shape-only - and never a
       value more favourable than the whole-member ratio (min of the two). */
    var vPtsKc = [];
    endsList().forEach(function (e) { if (e.uy) vPtsKc.push(+((e.x * 1000).toFixed(3))); });
    (S.ltbRestraints || []).forEach(function (r) { if (r.v !== false) { var xr = (+r.pos) * 1000; if (isFinite(xr) && xr >= -1e-6 && xr <= a.L + 1e-6) vPtsKc.push(+xr.toFixed(3)); } });
    vPtsKc = vPtsKc.filter(function (x, i) { return vPtsKc.indexOf(x) === i; }).sort(function (p, q) { return p - q; });
    var hasIntermediateKc = !isCant && vPtsKc.length >= 3;
    function bayC1For(ev) {
      if (!hasIntermediateKc || sec.isBox || !ev.sol.mode || !ev.sol.mode.x) return null;
      var m = ev.sol.mode, pkv = 0, px = null;
      m.x.forEach(function (x, i) { if (Math.abs(m.phi[i]) > pkv) { pkv = Math.abs(m.phi[i]); px = x; } });
      if (px == null) return null;
      var xa = vPtsKc[0], xb = vPtsKc[vPtsKc.length - 1];
      for (var i2 = 0; i2 < vPtsKc.length - 1; i2++) if (px >= vPtsKc[i2] - 1e-6 && px <= vPtsKc[i2 + 1] + 1e-6) { xa = vPtsKc[i2]; xb = vPtsKc[i2 + 1]; break; }
      var Ls = xb - xa; if (Ls < 1e-3) return null;
      var Mf = momentFromSamples(ev.res.fb.xs, ev.res.fb.M);
      var sp2 = secProps(sec);
      try {
        var cpl2 = unitLoadsFor(a, ev.res.combo).couples.filter(function (x) { return x > xa + 1e-9 && x < xb - 1e-9; }).map(function (x) { return x - xa; });
        var base2 = { E: a.E, G: G_STEEL, Iz: sp2.Iz, It: sp2.It, Iw: sp2.Iw, L: Ls, zj: 0,
                      restraints: [{ x: 0 }, { x: Ls }], nElem: 24, refine: true, couplePositions: cpl2 };
        var shb = mcrEigen(Object.assign({}, base2, { moment: (function (x0) { return function (x) { return Mf(x0 + x); }; })(xa) }));
        var unb = mcrEigen(Object.assign({}, base2, { moment: function () { return 1e6; } }));
        if (!(unb.Mcr > 0) || !shb.converged || !unb.converged) return null;
        return { C1: shb.Mcr / unb.Mcr, a: xa, b: xb };
      } catch (e2) { return null; }
    }

    /* kc = 1/sqrt(C1), NA 2.18. C1 is SHAPE-ONLY (zg = zj = 0): load height
       and monosymmetry are already inside Mcr, and folding them into kc
       would double-count them. Not applied to cantilevers (no published kc). */
    function chiChainFor(sol0, MxC, bay) {
      var lam = Math.sqrt(Wy * fy / sol0.Mcr);
      var C1c = sol0.C1;
      var c1lbl = 'back-calculated from the eigen solution: M<sub>cr</sub>(z<sub>g</sub>=0,z<sub>j</sub>=0)/M<sub>cr</sub>(uniform moment, same restraints)';
      if (bay && isFinite(bay.C1) && bay.C1 > 0 && bay.C1 < C1c) {
        C1c = bay.C1;
        c1lbl = 'critical bay ' + g(bay.a / 1000, 2) + '&ndash;' + g(bay.b / 1000, 2) + ' m between restraints (NA 2.18 basis): bay M<sub>cr</sub>(shape)/M<sub>cr</sub>(uniform); the whole-member ratio ' + sol0.C1.toFixed(3) + ' is not taken for k<sub>c</sub>';
      }
      var trusted = sol0.c1Converged;
      if (S.C1o != null) {
        c1lbl = 'user override for k<sub>c</sub> (eigen value was ' + C1c.toFixed(3) + ')';
        C1c = S.C1o;
        trusted = true;
      }
      /* k_c = 1/sqrt(C1) floored at 1/sqrt(2.76) = 0.60, the Table 6.6 lower
         bound (psi = -1); the back-calculated eigen C1 is unbounded (19 Sep
         2026 gap closure, item 3.5). kcFromC1() in eurocode-checks.js. */
      var kcc = 1.0, kcRaw = 1.0, kcFloored = false;
      if (trusted) {
        var kcr = kcFromC1(C1c);   // eurocode-checks.js is loaded before this patch (asserted above)
        kcc = kcr.kc; kcRaw = kcr.kcRaw; kcFloored = kcr.floored;
      }
      var Phi0 = null, chi0 = 1, f0 = 1, chiMod0 = 1, ign0 = true;
      if (lam > 0.4) {
        Phi0 = 0.5 * (1 + curve.alphaLT * (lam - 0.4) + 0.75 * lam * lam);
        chi0 = Math.min(1 / (Phi0 + Math.sqrt(Math.max(Phi0 * Phi0 - 0.75 * lam * lam, 1e-12))), 1, 1 / (lam * lam));
        f0 = isCant ? 1 : Math.min(1 - 0.5 * (1 - kcc) * (1 - 2 * Math.pow(lam - 0.8, 2)), 1);
        chiMod0 = Math.min(chi0 / f0, 1, 1 / (lam * lam));
        ign0 = false;
      }
      var MbRd0 = Math.min(chiMod0 * Wy * fy / gM1 / 1e6, b.McRd);
      return { lamLT: lam, C1: C1c, c1label: c1lbl, c1Trusted: trusted, kc: kcc, kcRaw: kcRaw, kcFloored: kcFloored,
               Phi: Phi0, chi: chi0, f: f0, chiMod: chiMod0, ign: ign0,
               MbRd: MbRd0, util: MbRd0 > 0 ? MxC / MbRd0 : 99, MxC: MxC };
    }

    var govEv = null, minMcrEv = null;
    evals.forEach(function (ev) {
      ev.bay = bayC1For(ev);
      ev.chain = chiChainFor(ev.sol, Math.abs(ev.res.Mmax) / 1e6, ev.bay);
      if (!govEv || ev.chain.util > govEv.chain.util) govEv = ev;
      if (!minMcrEv || ev.sol.Mcr < minMcrEv.sol.Mcr) minMcrEv = ev;
    });

    var sol = govEv.sol, chn = govEv.chain;
    var Mcr = sol.Mcr / 1e6;                     // kN.m
    var lamLT = chn.lamLT;
    var C1 = chn.C1, c1label = chn.c1label, c1Trusted = chn.c1Trusted, kc = chn.kc, kcRaw = chn.kcRaw, kcFloored = chn.kcFloored;
    if (kcFloored) warn.push('C<sub>1</sub> = ' + C1.toFixed(3) + ' exceeds 2.76: k<sub>c</sub> = 1/&radic;C<sub>1</sub> = ' + kcRaw.toFixed(3) + ' is floored at 1/&radic;2.76 = ' + kc.toFixed(3) + ', the Table 6.6 lower bound (&psi; = &minus;1 end-moment case), so the f-factor is not extrapolated beyond its calibrated range.');
    if (!c1Trusted) warn.push('The reference solves used to back-calculate C<sub>1</sub> did not converge; k<sub>c</sub> = 1.0 has been used, which is the conservative value (f = 1.0, hence the lower M<sub>b,Rd</sub>). M<sub>cr</sub> itself is unaffected.');

    var Phi = chn.Phi, chi = chn.chi, f = chn.f, chiMod = chn.chiMod, ign = chn.ign;
    var MbRd = chn.MbRd;

    /* Channel + eccentric load: the eigen Mcr is the LTB half of the story.
       The primary torque from e0 must be carried by the EN 1993-6 Annex A
       interaction. If that check cannot run, do not allow a PASS. */
    var chanTorsionGap = false;
    if (sec.kind === 'channel' && a.tors && a.tors.on && !(b.tor && b.tor.p385)) {
      chanTorsionGap = true;
      unsupported.push('Channel with eccentric load: the LTB check is valid (z<sub>j</sub> = 0, I<sub>w</sub> about the shear centre), ' +
        'but the primary torque from e<sub>0</sub> must be combined with it through the EN 1993-6 Annex A interaction, ' +
        'which is not available for this support/load arrangement (it needs a fork-fork single span with full-span or point torques). PASS is blocked.');
    }

    ltb = { eigen: true, na: false, cant: isCant, channel: sec.kind === 'channel', box: !!sec.isBox,
            Mcr: Mcr, McrRev: sol.McrRev / 1e6, McrUniform: sol.McrUniform / 1e6, McrShape: sol.McrShape / 1e6,
            C1: C1, c1label: c1label, kc: kc, kcRaw: kcRaw, kcFloored: kcFloored, lamLT: lamLT, lamLTmcr: lamLT,
            curve: curve, Phi: Phi, chi: chi, f: f, chiMod: chiMod, ign: ign, ignM: ign,
            chiM: chi, chiModM: chiMod, fM: f, PhiM: Phi,           // aliases: Annex A block reads chiM
             MbRd: MbRd, MbMcr: MbRd, MbSimp: MbRd, McrBack: Mcr,    // aliases: Annex A reads Mcr / McrBack
             zg: sol.zg, zj: sol.zj, nElem: sol.nElem, meshError: meshErrorMax,
             zgValues: sol.zgValues, zgUniform: sol.zgUniform,
             c1Trusted: c1Trusted, mcrConverged: allConverged,
             mode: sol.mode, warn: warn, chanTorsionGap: chanTorsionGap,
             MxGov: chn.MxC, governCombo: govEv.res.combo ? govEv.res.combo.label : '',
             nCombos: evals.length, nSolves: shared.stats.solved, nCached: shared.stats.cached,
             Iz: sol.sp.Iz, It: sol.sp.It, Iw: sol.sp.Iw, hs: sol.sp.hs };

    /* ---- Informational extras (NOT the design basis) ----
       1. Locate the buckled-mode twist peak: shows which bay the eigenmode
          localises in - the "critical segment" of the conventional method.
       2. Conventional segment-method comparison: each bay between lateral
          (v) restraint points solved in ISOLATION with fork ends and its own
          share of the moment diagram/loads. This discards lateral-bending and
          warping continuity across the restraints, so it is conservative
          relative to the whole-member eigen Mcr, which remains the design
          basis. Computed only when intermediate lateral restraints exist. */
    var modePeakX = null;
    if (sol.mode && sol.mode.x && sol.mode.x.length) {
      var mpk = 0;
      sol.mode.x.forEach(function (x, i) { if (Math.abs(sol.mode.phi[i]) > mpk) { mpk = Math.abs(sol.mode.phi[i]); modePeakX = x; } });
    }
    ltb.modePeakX = modePeakX;
    var vPts = [];
    endsList().forEach(function (e) { if (e.uy) vPts.push(+((e.x * 1000).toFixed(3))); });
    (S.ltbRestraints || []).forEach(function (r) { if (r.v !== false) { var xr = (+r.pos) * 1000; if (isFinite(xr) && xr >= -1e-6 && xr <= a.L + 1e-6) vPts.push(+xr.toFixed(3)); } });
    vPts = vPts.filter(function (x, i) { return vPts.indexOf(x) === i; }).sort(function (p, q) { return p - q; });
    ltb.vPoints = vPts;
    ltb.segments = null;
    if (!isCant && vPts.length >= 3) {
      var gfbG = govEv.res.fb, MfunG = momentFromSamples(gfbG.xs, gfbG.M);
      var ulG = unitLoadsFor(a, govEv.res.combo);
      var spG = secProps(sec), zjG = zjFor(sec);
      var segs = [];
      for (var si = 0; si < vPts.length - 1; si++) {
        var xa = vPts[si], xb = vPts[si + 1], Ls = xb - xa;
        if (Ls < 1e-3) continue;
        var dlS = ulG.distLoads.map(function (d) {
          var x1 = Math.max(d.x1, xa), x2 = Math.min(d.x2, xb);
          if (x2 - x1 < 1e-9) return null;
          var wAt = function (x) { return (d.x2 - d.x1 < 1e-9) ? d.w1 : d.w1 + (d.w2 - d.w1) * (x - d.x1) / (d.x2 - d.x1); };
          return { x1: x1 - xa, x2: x2 - xa, w1: wAt(x1), w2: wAt(x2), zg: d.zg || 0 };
        }).filter(function (d) { return d; });
        var plS = ulG.pointLoads.filter(function (q) { return q.x > xa + 1e-9 && q.x < xb - 1e-9; })
          .map(function (q) { return { x: q.x - xa, P: q.P, zg: q.zg || 0 }; });
        var cplS = ulG.couples.filter(function (x) { return x > xa + 1e-9 && x < xb - 1e-9; }).map(function (x) { return x - xa; });
        var seg = { a: xa, b: xb, ok: false };
        try {
          var rsS = mcrEigen({ E: a.E, G: G_STEEL, Iz: spG.Iz, It: spG.It, Iw: spG.Iw, L: Ls, zj: zjG,
            restraints: [{ x: 0 }, { x: Ls }], moment: (function (x0) { return function (x) { return MfunG(x0 + x); }; })(xa),
            distLoads: dlS, pointLoads: plS, couplePositions: cplS, nElem: 32, refine: true });
          var MsS = 0; for (var ii = 0; ii <= 200; ii++) { var mmS = Math.abs(MfunG(xa + Ls * ii / 200)); if (mmS > MsS) MsS = mmS; }
          var lamSg = Math.sqrt(Wy * fy / rsS.Mcr);
          var chiSg = 1;
          if (lamSg > 0.4) {
            var PhiSg = 0.5 * (1 + curve.alphaLT * (lamSg - 0.4) + 0.75 * lamSg * lamSg);
            chiSg = Math.min(1 / (PhiSg + Math.sqrt(Math.max(PhiSg * PhiSg - 0.75 * lamSg * lamSg, 1e-12))), 1, 1 / (lamSg * lamSg));
          }
          var MbSg = Math.min(chiSg * Wy * fy / gM1 / 1e6, b.McRd);
          seg = { a: xa, b: xb, ok: true, Ms: MsS / 1e6, Mcr: rsS.Mcr / 1e6, lam: lamSg, chi: chiSg, Mb: MbSg,
                  util: MbSg > 0 ? (MsS / 1e6) / MbSg : 99 };
        } catch (eS) { seg.err = eS.message; }
        segs.push(seg);
      }
      if (segs.length) ltb.segments = segs;
    }

    var Mx = b.Mx;
    var ltbUtil = chn.util;   // governing combination's own Mmax / its MbRd
    var ltbBasis = 'elastic critical moment from the finite-element eigenvalue solution ' +
      '(4 DOF/node: v, v\', &phi;, &phi;\'; z<sub>g</sub> and z<sub>j</sub> included; ' + ltb.nElem + ' elements, Richardson-extrapolated' +
      (evals.length > 1 ? '; ' + evals.length + ' ULS combinations each solved with their own moment diagram, governing: ' + ltb.governCombo : '') + ')';

    /* ---- BAY ISOLATION between intermediate lateral restraints ----
       With intermediate lateral restraints (>= 2 bays between points holding
       v) each bay is also checked in ISOLATION with fork ends, its own Mcr and
       its own peak moment - the conventional segment method - and the LTB
       utilisation is the MAX of the worst isolated bay and the whole-member
       eigen result, so a bay is never under-checked; isolating a bay discards
       the v' / phi' continuity across the restraint and is conservative
       relative to the whole-member Mcr. ltb.segments is empty without
       intermediate restraints, so an unrestrained span is unaffected. */
    var spanGov = null;
    if (ltb.segments && ltb.segments.length) {
      ltb.segments.forEach(function (s2) { if (s2.ok && (!spanGov || s2.util > spanGov.util)) spanGov = s2; });
    }
    if (spanGov && spanGov.util >= ltbUtil) {
      ltbUtil = spanGov.util;
      ltb.spanGoverns = true;
      ltb.spanGov = { a: spanGov.a, b: spanGov.b, Ms: spanGov.Ms, Mcr: spanGov.Mcr, lam: spanGov.lam, chi: spanGov.chi, Mb: spanGov.Mb, util: spanGov.util };
      ltbBasis = 'intermediately restrained member: lateral-torsional buckling checked BAY BY BAY between the lateral restraints. ' +
        'Each bay is solved in isolation with fork ends, its own M<sub>cr</sub> and its own peak moment; the worst bay governs. ' +
        'Governing bay ' + g(spanGov.a / 1000, 2) + '&ndash;' + g(spanGov.b / 1000, 2) + ' m: M<sub>Ed</sub> = ' + f1(spanGov.Ms, 1) + ', M<sub>cr</sub> = ' + f1(spanGov.Mcr, 1) + ', M<sub>b,Rd</sub> = ' + f1(spanGov.Mb, 1) + ' kN&middot;m, utilisation ' + g(spanGov.util, 2) + '. ' +
        '(Whole-member eigen M<sub>cr</sub> = ' + f1(ltb.Mcr, 1) + ' kN&middot;m retained for reference.)';
    }

    /* ---- Standard-method comparison (no eigen solve): closed-form Mcr for
       the same segment - the governing bay when the bay-by-bay check
       governs, otherwise the whole member - and for the SAME COMBINATION as
       the design eigen value (govEv: its own moment diagram for C1 and the
       reference software-style C1 inputs, its own load factors for the SN003a shape
       recognition and the load height z_g = most destabilising per-load value
       of that combination), with C1 from sn003aC1 (SN006a for a cantilever),
       LE = LE-factor x segment length and C2 where published. ---- */
    var segC1 = spanGov ? { xa: spanGov.a, xb: spanGov.b, whole: false }
                        : (typeof c1Segment === 'function' ? c1Segment(a) : { xa: 0, xb: a.L, whole: true });
    var stdCmp = null;
    try { stdCmp = mcrStandardFor(a, sec, segC1, { fb: govEv.res.fb, factors: govEv.res.combo.factors, combo: govEv.res.combo }); }
    catch (eStd) { stdCmp = { route: 'n/a', Mcr: null, C1: null, label: 'closed form not available: ' + eStd.message, c1in: null, seg: segC1 }; }
    if (stdCmp && !stdCmp.c1in && typeof c1Inputs === 'function') stdCmp.c1in = c1Inputs(govEv.res.fb, segC1.xa, segC1.xb);
    ltb.mcrMethod = 'eigen';
    ltb.McrEigen = spanGov ? spanGov.Mcr : Mcr;          // kN.m, the design value's segment
    ltb.McrStandard = (stdCmp && stdCmp.Mcr != null && isFinite(stdCmp.Mcr)) ? stdCmp.Mcr : null;
    ltb.McrRatio = (ltb.McrStandard > 0) ? ltb.McrEigen / ltb.McrStandard : null;
    ltb.std = stdCmp;
    ltb.c1in = stdCmp ? stdCmp.c1in : null;
    ltb.c1seg = segC1;
    ltb.c1route = stdCmp ? stdCmp.route : null;

    /* ---- EN 1993-6 Annex A: LTB + minor-axis bending + torsion ---- */
    var annex = null;
    if (b.tor && b.tor.p385) {
      /* The Annex A amplifier k_alpha uses the SMALLEST Mcr across the
         evaluated ULS combinations (conservative when several diagram shapes
         exist), with the matching chi_LT (no f-factor - P385 basis). */
      var lamAA = Math.sqrt(Wy * fy / minMcrEv.sol.Mcr);
      var chiA = 1;
      if (lamAA > 0.4) {
        var PhiAA = 0.5 * (1 + curve.alphaLT * (lamAA - 0.4) + 0.75 * lamAA * lamAA);
        chiA = Math.min(1 / (PhiAA + Math.sqrt(Math.max(PhiAA * PhiAA - 0.75 * lamAA * lamAA, 1e-12))), 1, 1 / (lamAA * lamAA));
      }
      var MbA = chiA * Wy * fy / gM1 / 1e6;
      var McrA = minMcrEv.sol.Mcr / 1e6;
      /* Cmz is an equivalent-uniform-moment factor for the minor-axis moment
         diagram (EN 1993-1-1 Annex B, Table B.3). C1 describes the major-axis
         diagram's effect on Mcr. Different quantities, different diagrams.

         The old lookup only worked because the retired sn003aC1() returned C1
         from a discrete table, where 1.348 meant "SS + central point load" and
         1.127 meant "SS + UDL". It was a proxy for the load case, not for C1.
         C1 is now a continuous eigenvalue ratio that also absorbs intermediate
         restraints and the end conditions, so hitting a +/-0.02 window is
         coincidental and can reduce the minor-axis demand spuriously.

         1.0 is conservative: the term enters additively as +Cmz*Mz/MzR.
         Proper derivation = Table B.3 applied to the Mz(x)=phi(x)*My(x)
         diagram in b.tor.grids. NOTE this is not cmTableB3(), which reads the
         major-axis diagram a.governM.fb and returns Cmy. Not implemented.
         20 Sep 2026 review: the override S.Cmzo applies to the twist-induced
         diagram only - with an imposed M_z (a constant diagram, psi = 1) Table
         B.3 gives C_mz = 1.0 and annexAEval() uses 1.0 whatever the override
         (annex.CmzBasis says so; no UI input exists for S.Cmzo - a state
         field only, js/03-state-ui.js). */
      var Cmz = (S.Cmzo != null && isFinite(S.Cmzo)) ? +S.Cmzo : 1.0;
      /* 20 Sep 2026 torsion + N/Mz: Eq (A.1) per station through the shared
         annexAEval() of js/checks/eurocode-checks.js - M_z,Ed = M_z,tot =
         imposed M_z + phi.M_y in both the C_mz M_z/M_z,Rd term and k_zw;
         C_mz = 1.0 whenever an M_z is imposed (the override applies to the
         twist-induced diagram only); identical to the former inline loop
         when no M_z is imposed. */
      annex = annexAEval(b.tor, MbA, McrA, Cmz);
    }

    var useB1u = sec.isBox || (ltb.MbRd >= b.McRd * 0.9999);
    // Bound the interaction with the lowest resistance of all enabled diagrams.
    // annexB2 separately sweeps each combination's own moment and Cm.
    var memberMb = Math.min.apply(null, evals.map(function(ev){ return ev.chain.MbRd; }));
    if (spanGov) memberMb = Math.min(memberMb, spanGov.Mb);
    useB1u = sec.isBox || (memberMb >= b.McRd * 0.9999);
    var buck = (b.ax && !b.ax.tension) ? annexB2(a, sec, fy, b.cl, memberMb, useB1u, isCant, b.aeff) : null;
    if (buck && buck.tfb && !buck.tfb.ok) unsupported.push('PFC under axial compression: ' + buck.tfb.reason + '; torsional / torsional-flexural buckling (cl 6.3.1.4) cannot be verified, PASS is blocked.');
    var restraintF = restraintForces(a, sec);
    if (buck && buck.lczFromRestraints)
      warn.push('Minor-axis strut buckling length L<sub>cr,z</sub> = ' + (buck.LcrZ / 1000).toFixed(2) + ' m, taken as the largest spacing between adjacent lateral restraint points (SCI P360 6.2: secondary members act as bracing points; k = 1.0 between restraints). ' +
        'Ensure each restraint really is an effective bracing point - adequate stiffness, strength and anchorage. The major axis keeps L<sub>cr,y</sub> = L<sub>E</sub>&times;L = ' + (buck.LcrY / 1000).toFixed(2) + ' m.');

    var utils = [
      { name: 'Shear  V_Ed/V_c,Rd', val: b.shearUtil },
      { name: 'Bending  M_Ed/M_c,Rd', val: b.momUtil },
      { name: 'LTB  M_Ed/M_b,Rd', val: ltbUtil },
      { name: 'Deflection', val: b.dmax / b.dlimit }
    ];
    if (b.ax) {
      utils.push({ name: b.ax.tension ? 'Tension  N_Ed/N_t,Rd' : ((b.ax.aeff && b.ax.aeff.active) ? 'Compression  N_Ed/N_c,Rd (A_eff)' : 'Compression  N_Ed/N_pl,Rd'), val: b.ax.nUtil });
      utils.push({ name: b.ax.biax ? ('Biaxial bending' + ((S.axial || 0) !== 0 ? ' + axial' : '') + ' (6.2.9.1)') : 'Bending+axial cross-section (6.2.9)', val: b.ax.mUtil });
      /* Eq 6.61/6.62 are needed with axial compression AND for biaxial bending
         on an LTB-susceptible member with N_Ed = 0: Eq 6.62 then reads
         kzy*My/Mb,Rd + kzz*Mz/Mcz,Rd, which the separate LTB and cross-section
         checks do not cover. */
      if (!b.ax.tension && buck && (buck.Fc > 1e-9 || buck.biax)) {
        if (buck.tfb && buck.tfb.ok) utils.push({ name: TFB_UTIL_NAME, val: buck.tfb.util });
        utils.push({ name: 'Member buckling y-y (Eq 6.61)', val: buck.u1 });
        utils.push({ name: 'Member buckling z-z (Eq 6.62)', val: buck.u2 });
      }
    }
    if (annex) utils.push({ name: 'LTB+torsion (EN 1993-6 Annex A)', val: annex.u });
    if (b.coex) utils.push({ name: b.coex.pureShearFail ? 'Pure shear failure at M-V check point (6.2.6)' : 'Bending+shear coexistent (6.2.8)', val: b.coex.u });
    if (b.mvn) utils.push({ name: mvnUtilName(b.mvn), val: b.mvn.u });
    if (b.web && b.web.checked) {
      utils.push({ name: WEB_UTIL_NAMES[0], val: b.web.util2 });
      utils.push({ name: WEB_UTIL_NAMES[1], val: b.web.util72 });
    }
    if (b.tor && b.tor.box) {
      utils.push({ name: 'Torsion  T_Ed/T_Rd', val: b.tor.torUtil });
      utils.push({ name: 'Shear+torsion  V_Ed/V_pl,T,Rd', val: b.tor.vtUtil });
    }
    if (b.tor && b.tor.p385) {
      utils.push({ name: 'Bending+torsion cross-section (P385 3.1.2)', val: b.tor.cross.u });
      utils.push({ name: 'Shear+torsion  V_Ed/V_pl,T,Rd', val: b.tor.vtUtil });
    }
    /* 20 Sep 2026 torsion + N/Mz: elastic yield criterion (6.1), cl 6.2.7(5)
       (ELASTIC_TORSION_UTIL_NAME, eurocode-checks.js); 20 Sep 2026 review: in utils
       only when verdict-binding (elasticBindingPolicy of checksEC3Restrained: Class 3,
       or an open Class 1/2 section with N_Ed), otherwise in b.info (carried over from
       the restrained result with its advisory when it exceeds 1) */
    if (b.tor && b.tor.elastic && b.tor.elastic.binding) utils.push({ name: ELASTIC_TORSION_UTIL_NAME, val: b.tor.elastic.u });
    /* 20 Sep 2026 torsion + N/Mz: advisory superposition of Eq 6.62 and (A.1) with
       this path's k_alpha (information only, never a utilisation) */
    var advisory = (b.advisory || []).slice();
    if (annex && annex.unbounded) advisory.push(KALPHA_UNBOUNDED_FAIL);   // 20 Sep 2026: M_y,Ed >= M_cr is a FAIL row (utilisation 99), not NOT VERIFIED
    if (b.tor && b.tor.p385) { b.tor.superposition = torsionSuperposition(b.tor, annex, buck); if (b.tor.superposition) advisory.push(b.tor.superposition.text); }
    /* 20 Sep 2026 review: the basis text of this path (torsionCombinedBasis, eurocode-checks.js) */
    if (b.tor && (b.tor.p385 || b.tor.box)) b.tor.combinedBasis = torsionCombinedBasis({ box: !!b.tor.box, tension: !!(b.ax && b.ax.tension), hasN: Math.abs(S.axial || 0) > 1e-9, hasMz: Math.abs(S.Mz || 0) > 1e-9, restrained: false, buckEvaluated: !!(!(b.ax && b.ax.tension) && buck && (buck.Fc > 1e-9 || buck.biax)), annexEvaluated: !!annex, binding: !!(b.tor.elastic && b.tor.elastic.binding) });
    var gov = utils[0]; utils.forEach(function (u) { if (u.val > gov.val) gov = u; });
    var pass = unsupported.length === 0 && utils.every(function (u) { return u.val <= 1.0001; });

    return Object.assign({}, b, { sci: false, sciU: true, mcrMethod: 'eigen', unsupported: unsupported, advisory: advisory, ltb: ltb,
      ltbUtil: ltbUtil, ltbBasis: ltbBasis, C1: C1, c1label: c1label, LE: a.L,
      utils: utils, gov: gov, pass: pass, annex: annex, buck: buck, restraintForces: restraintF });
  };

  /* ================================================================
     PART 5 - report block
     Replaces the `sciUltbBlocks` template in render().
     Pure ASCII + HTML entities, so it survives the file's cp1252 encoding.
     ================================================================ */
  window.ltbEigenReport = function (c, a, sec) {
    var LT = c.ltb || {};
    if (LT.na && LT.closed) return '<div class="section-title smallgap">Lateral&ndash;Torsional Buckling (Cl. 6.3.2.1(2))</div>' +
      '<div class="calc-block">' +
      '<div>Closed hollow section</div><div class="formula">' + (S.family === 'rhs' ? 'RHS' : 'SHS') + ' / closed box section &mdash; not susceptible to lateral-torsional buckling</div><div class="value">LTB not required</div><div class="status ok">Not required</div>' +
      '<div>M<sub>b,Rd</sub> = M<sub>c,Rd</sub></div><div class="formula">Full cross-section bending resistance used directly; adequacy is governed by the Clause 6.2 moment check above</div><div class="value">' + f1(LT.MbRd, 1) + ' kN&middot;m</div><div></div>' +
      '</div>';
    if (!LT.eigen) return '';
    if (LT.failed) return '<div class="section-title smallgap">Lateral&ndash;Torsional Buckling</div>' +
      '<div class="calc-block"><div>M<sub>cr</sub></div><div class="formula">' + LT.err +
      '</div><div class="value">&mdash;</div><div class="status fail">BLOCKED</div></div>';

    var Wy = c.cl.cls <= 2 ? sec.Sx : sec.Zx;
    var restr = (S.ltbRestraints || []).length;
    var rows = '';
    rows += '<div>Buckling model</div><div class="formula">FE eigenvalue solution of (K<sub>e</sub> + &lambda;K<sub>g</sub>)d = 0 over the governing BMD; ' +
            'Hermite cubics, 4 DOF/node (v, v&prime;, &phi;, &phi;&prime;); ' + LT.nElem + ' elements, Richardson-extrapolated</div>' +
            '<div class="value">mesh err &lt; ' + g(Math.max(LT.meshError || 0, 1e-6) * 100, 3) + ' %</div><div></div>';
    rows += '<div>Section properties</div><div class="formula">I<sub>z</sub> = ' + g(LT.Iz / 1e4, 0) + ' cm<sup>4</sup>; I<sub>T</sub> = ' + g(LT.It / 1e4, 1) +
            ' cm<sup>4</sup>; I<sub>w</sub> = ' + g(LT.Iw / 1e12, 4) + ' dm<sup>6</sup>' + (LT.channel ? ' (about the shear centre)' : '') +
            '; G = 81000 N/mm&sup2;</div><div class="value">z<sub>j</sub> = ' + g(LT.zj, 1) + ' mm</div><div></div>';
    var endBc = ltbEndBcText(S, sec);
    rows += '<div>End boundary conditions</div><div class="formula">' + endBc +
            (restr ? '; ' + restr + ' intermediate restraint(s)' : '') +
            '</div><div class="value">&mdash;</div><div></div>';
    var zref = (typeof loadHeightReference === 'function') ? loadHeightReference(sec) : null;
    var zrefText = zref ? '; refs: top +' + g(zref.topSurface, 0) + ' mm, bottom ' + g(zref.bottomSurface, 0) + ' mm' : '';
    var zfmt = function(z){ return (z>0?'+':'') + g(z,0); };
    var zgText = (LT.zgValues && LT.zgValues.length > 1)
      ? 'per-load z<sub>g</sub> = ' + LT.zgValues.map(zfmt).join(', ') + ' mm above the shear centre'
      : 'z<sub>g</sub> = ' + zfmt(LT.zg || 0) + ' mm above the shear centre';
    rows += '<div>Load height</div><div class="formula">' + zgText + zrefText +
            (LT.zg > 0 ? ' (max value destabilising)' : LT.zg < 0 ? ' (max value stabilising)' : '') + '</div><div class="value">M<sub>cr</sub> (load reversed) = ' + f1(LT.McrRev, 1) + ' kN&middot;m</div><div></div>';
    rows += '<div><b>M<sub>cr</sub></b> (FE eigenvalue method)</div><div class="formula">eigenvalue &times; max|M(x)| &mdash; no C<sub>1</sub>, C<sub>2</sub> or C<sub>3</sub> used' +
            (LT.nCombos > 1 ? '; each of the ' + LT.nCombos + ' ULS combinations solved with its own diagram &mdash; governing: ' + LT.governCombo : '') +
            (LT.nSolves != null ? '; ' + LT.nSolves + ' eigen solve(s) this render' + (LT.nCached ? ', ' + LT.nCached + ' from the cache' : '') : '') + '</div>' +
            '<div class="value"><b>' + f1(LT.Mcr, 1) + ' kN&middot;m</b></div><div></div>';
    /* standard closed-form comparison for the same segment (informational) */
    if (LT.std) {
      var ci = LT.c1in, sg = LT.c1seg || {};
      var segTxt = (sg.whole === false) ? 'segment ' + g(sg.xa / 1000, 2) + '&ndash;' + g(sg.xb / 1000, 2) + ' m' : 'whole member';
      var ciTxt = ci ? 'M<sub>1</sub> = ' + f1(ci.M1, 1) + ', M<sub>2</sub> = ' + f1(ci.M2, 1) + ', M<sub>o</sub> = ' + f1(ci.Mo, 1) + ' kN&middot;m; &psi; = ' + f1(ci.psi, 3) + '; &mu; = ' + f1(ci.mu, 3) + '; ' : '';
      rows += '<div>Standard closed-form M<sub>cr</sub> (comparison, not the design basis)</div><div class="formula">' +
              (LT.std.route === 'sn006a' ? 'SN006a cantilever: C = ' : 'SN003a: C<sub>1</sub> = fn(M<sub>1</sub>, M<sub>2</sub>, M<sub>o</sub>, &psi;, &mu;) = ') +
              (LT.std.C1 != null ? g(LT.std.C1, 3) : '&mdash;') + ' &mdash; ' + ciTxt + (LT.std.label || '') +
              (LT.std.LE ? '; L<sub>E</sub> = ' + g(LT.std.LE / 1000, 2) + ' m (' + segTxt + ')' : '') +
              (LT.std.zgUsed ? '; C<sub>2</sub>z<sub>g</sub> term applied' : '') + '</div>' +
              '<div class="value">' + (LT.McrStandard != null ? f1(LT.McrStandard, 1) + ' kN&middot;m' : 'not covered') +
              (LT.McrRatio != null ? '<br>eigen / standard = ' + f1(LT.McrRatio, 2) : '') + '</div><div></div>';
    }
    rows += '<div>&lambda;&#772;<sub>LT</sub> = &radic;(W<sub>y</sub>f<sub>y</sub>/M<sub>cr</sub>)</div><div class="formula">&radic;(' + g(Wy, 0) + '&times;10&sup3;&times;' + g(a.py, 0) + '/' + f1(LT.Mcr, 1) + '&times;10<sup>6</sup>)</div><div class="value">' + f1(LT.lamLT, 3) + '</div><div></div>';
    rows += '<div>Buckling curve</div><div class="formula">' + (sec.isBox ? 'closed section, not listed in NA Table 6.3' : sec.kind === 'channel' ? 'not doubly symmetric' : 'NA Table 6.3, h/b = ' + g(sec.D / sec.B, 2)) +
            '</div><div class="value">curve ' + LT.curve.curve + ' (&alpha;<sub>LT</sub> = ' + g(LT.curve.alphaLT, 2) + ')</div><div></div>';
    if (LT.ign) {
      rows += '<div>&lambda;&#772;<sub>LT</sub> &le; 0.4 (NA 2.17)</div><div class="formula">LTB effects may be ignored (cl 6.3.2.2(4))</div><div class="value">&chi;<sub>LT,mod</sub> = 1.000</div><div class="status ok">Ignored</div>';
    } else {
      rows += '<div>&Phi;<sub>LT</sub>; &chi;<sub>LT</sub></div><div class="formula">&lambda;&#772;<sub>LT,0</sub> = 0.4, &beta; = 0.75 (NA 2.17); &Phi; = ' + g(LT.Phi, 3) + '</div><div class="value">&chi;<sub>LT</sub> = ' + g(LT.chi, 3) + '</div><div></div>';
      rows += '<div>C<sub>1</sub> (for k<sub>c</sub> only)</div><div class="formula">' + LT.c1label + '</div><div class="value">C<sub>1</sub> = ' + g(LT.C1, 3) + '</div><div></div>';
      if (LT.cant) rows += '<div>k<sub>c</sub> / f</div><div class="formula">not applied to cantilevers (no published k<sub>c</sub>)</div><div class="value">f = 1.000</div><div></div>';
      else rows += '<div>k<sub>c</sub> = 1/&radic;C<sub>1</sub> &ge; 1/&radic;2.76; f = 1&minus;0.5(1&minus;k<sub>c</sub>)[1&minus;2(&lambda;&#772;<sub>LT</sub>&minus;0.8)&sup2;] &le; 1</div><div class="formula">' + (LT.kcFloored ? '1/&radic;' + g(LT.C1, 3) + ' = ' + g(LT.kcRaw, 3) + ' &rarr; floored at 1/&radic;2.76 = ' + g(LT.kc, 3) + ' (Table 6.6 lower bound, &psi; = &minus;1)' : 'k<sub>c</sub> = ' + g(LT.kc, 3) + ' (NA 2.18; floor 0.60 not reached)') + '</div><div class="value">f = ' + g(LT.f, 3) + '</div><div></div>';
      rows += '<div>&chi;<sub>LT,mod</sub> = &chi;<sub>LT</sub>/f &le; min(1, 1/&lambda;&#772;&sup2;)</div><div class="formula">' + g(LT.chi, 3) + '/' + g(LT.f, 3) + '</div><div class="value">' + g(LT.chiMod, 3) + '</div><div></div>';
    }
    rows += '<div>M<sub>b,Rd</sub> = &chi;<sub>LT,mod</sub>W<sub>' + (c.cl.cls <= 2 ? 'pl' : 'el') + ',y</sub>f<sub>y</sub>/&gamma;<sub>M1</sub> &le; M<sub>c,Rd</sub></div>' +
            '<div class="formula">' + g(LT.ign ? 1 : LT.chiMod, 3) + '&times;' + g(Wy, 0) + '&times;' + g(a.py, 0) + '/1.0</div><div class="value">' + f1(LT.MbRd, 1) + ' kN&middot;m</div><div></div>';
    var MxLTB = LT.spanGoverns ? LT.spanGov.Ms : ((LT.MxGov != null) ? LT.MxGov : c.Mx);
    var MbLTB = LT.spanGoverns ? LT.spanGov.Mb : LT.MbRd;
    var mEdLbl = LT.spanGoverns ? ' (governing bay ' + g(LT.spanGov.a / 1000, 2) + '&ndash;' + g(LT.spanGov.b / 1000, 2) + ' m)' : (LT.nCombos > 1 ? ' (governing combination)' : '');
    rows += '<div>M<sub>Ed</sub> / M<sub>b,Rd</sub>' + mEdLbl + '</div><div class="formula">' + f1(MxLTB, 1) + ' / ' + f1(MbLTB, 1) + (LT.spanGoverns ? '  &mdash; per bay (between lateral restraints)' : '') + '</div><div class="value">' + g(MxLTB / Math.max(MbLTB, 1e-9), 2) + '</div>' + st(c.ltbUtil <= 1, 'OK', 'exceeded');

    var warnHtml = (LT.warn && LT.warn.length) ? '<div class="note" style="margin-left:0">' + LT.warn.map(function (w) { return '&bull; ' + w; }).join('<br>') + '</div>' : '';

    /* ---- buckled mode shape (normalised) ---- */
    var modeHtml = '';
    var md = LT.mode;
    if (md && md.x && md.x.length > 2) {
      var Wp = 640, Hp = 150, pdd = 12, Lm = md.x[md.x.length - 1];
      var vmx = 0; md.v.forEach(function (vv) { vmx = Math.max(vmx, Math.abs(vv)); });
      var Xf = function (x) { return pdd + (Wp - 2 * pdd) * x / Lm; };
      var Yc = Hp / 2, ampl = Hp / 2 - pdd - 16;
      var pphi = '', pv = '';
      md.x.forEach(function (x, i) {
        pphi += (i ? ' L ' : 'M ') + Xf(x).toFixed(1) + ' ' + (Yc - ampl * md.phi[i]).toFixed(1);
        pv += (i ? ' L ' : 'M ') + Xf(x).toFixed(1) + ' ' + (Yc - ampl * (vmx > 0 ? md.v[i] / vmx : 0)).toFixed(1);
      });
      var marks = (LT.vPoints || []).map(function (x) {
        var xx = Xf(x);
        return '<line x1="' + xx.toFixed(1) + '" y1="' + pdd + '" x2="' + xx.toFixed(1) + '" y2="' + (Hp - pdd) + '" stroke="#b91c1c" stroke-width="1" stroke-dasharray="3,3"/>' +
               '<text x="' + xx.toFixed(1) + '" y="' + (Hp - 2) + '" font-size="9" text-anchor="middle" fill="#b91c1c">' + g(x / 1000, 2) + '</text>';
      }).join('');
      modeHtml =
        '<div class="section-title smallgap">Buckled Mode Shape &mdash; critical eigenmode (normalised)</div>' +
        '<div style="border:1px solid #d1d5db;border-radius:6px;padding:6px;background:#fff">' +
        '<svg viewBox="0 0 ' + Wp + ' ' + Hp + '" style="width:100%;max-width:' + Wp + 'px;display:block" xmlns="http://www.w3.org/2000/svg">' +
        '<line x1="' + pdd + '" y1="' + Yc + '" x2="' + (Wp - pdd) + '" y2="' + Yc + '" stroke="#9ca3af" stroke-width="1"/>' +
        marks +
        '<path d="' + pphi + '" fill="none" stroke="#1d4ed8" stroke-width="2"/>' +
        '<path d="' + pv + '" fill="none" stroke="#059669" stroke-width="1.6" stroke-dasharray="6,4"/>' +
        '</svg>' +
        '<div class="note" style="margin-left:0">Blue solid: twist &phi;(x); green dashed: lateral displacement v(x), each normalised to its own peak. ' +
        'Red dashed verticals: lateral restraint points. Peak twist at x = ' + (LT.modePeakX != null ? g(LT.modePeakX / 1000, 2) : '&mdash;') + ' m ' +
        '&mdash; the eigenmode localises in the critical bay, so the whole-member M<sub>cr</sub> already embodies the worst segment.</div></div>';
    }

    /* ---- conventional segment-method comparison (informational) ---- */
    var segHtml = '';
    if (LT.segments && LT.segments.length) {
      var worstSeg = null; LT.segments.forEach(function (s2) { if (s2.ok && (!worstSeg || s2.util > worstSeg.util)) worstSeg = s2; });
      var isBasis = !!LT.spanGoverns;
      segHtml =
        '<div class="section-title smallgap">' + (isBasis
          ? 'Lateral&ndash;Torsional Buckling &mdash; Bay by Bay (design basis, between lateral restraints)'
          : 'Segment-Method Comparison (informational &mdash; NOT the design basis)') + '</div>' +
        '<div class="calc-block">' +
        LT.segments.map(function (s2) {
          if (!s2.ok) return '<div>Bay ' + g(s2.a / 1000, 2) + '&ndash;' + g(s2.b / 1000, 2) + ' m</div><div class="formula">could not be solved in isolation: ' + (s2.err || '') + '</div><div class="value">&mdash;</div><div></div>';
          var isWorst = (worstSeg && s2 === worstSeg);
          var tag = isWorst ? ('<div class="status' + (isBasis ? (s2.util <= 1 ? ' ok' : ' fail') : '') + '">' + (isBasis ? 'governs' : 'worst segment') + '</div>') : '<div></div>';
          return '<div>Bay ' + g(s2.a / 1000, 2) + '&ndash;' + g(s2.b / 1000, 2) + ' m (between lateral restraints, fork ends)</div>' +
            '<div class="formula">M<sub>Ed</sub> = ' + f1(s2.Ms, 1) + '; M<sub>cr</sub> = ' + f1(s2.Mcr, 0) + '; M<sub>b,Rd</sub> = ' + f1(s2.Mb, 0) + ' kN&middot;m; &lambda;&#772;<sub>LT</sub> = ' + g(s2.lam, 3) + '; &chi;<sub>LT</sub> = ' + g(s2.chi, 3) + ' (no f-factor)</div>' +
            '<div class="value">util = ' + g(s2.util, 3) + '</div>' + tag;
        }).join('') +
        '</div>' +
        '<div class="note" style="margin-left:0">' + (isBasis
          ? 'Intermediately restrained member: LTB is checked <b>bay by bay</b>. Each bay between lateral restraints is isolated with fork ends and its own share of the moment diagram and loads, giving its own M<sub>cr</sub>; the worst bay governs the LTB utilisation above. Isolation conservatively discards the lateral-bending (v&prime;) and warping (&phi;&prime;) continuity across the restraints (whole-member eigen M<sub>cr</sub> = ' + f1(LT.Mcr, 1) + ' kN&middot;m is shown for reference only).'
          : 'Conventional check for comparison: each bay between lateral restraints is isolated with fork ends and its own share of the moment diagram and loads. ' +
            'Isolation discards the lateral-bending (v&prime;) and warping (&phi;&prime;) continuity that the adjacent bays provide across the restraint points, so these values are conservative. ' +
            'The design basis is the whole-member eigen solution above: M<sub>cr</sub> = ' + f1(LT.Mcr, 1) + ' kN&middot;m' + (worstSeg ? ' vs the worst isolated segment ' + f1(worstSeg.Mcr, 0) + ' kN&middot;m' : '') + '.') + '</div>';
    }

    return '<div class="section-title smallgap">Lateral&ndash;Torsional Buckling &mdash; Elastic Critical Moment (FE eigenvalue solution)</div>' +
           '<div class="calc-block">' + rows + '</div>' + warnHtml + modeHtml + segHtml;
  };

  /* ================================================================
     PART 6 - intermediate lateral restraint UI
     ================================================================ */
  if (S.ltbRestraints == null) S.ltbRestraints = [];
  if (S.zj == null) S.zj = 0;
  if (S.Cmzo === undefined) S.Cmzo = null;   // verified Cmz override; null = conservative 1.0

  function injectUI() {
    var host = document.getElementById('restraintRow');
    if (!host || document.getElementById('ltbRestraintPanel')) return;
    var div = document.createElement('div');
    div.id = 'ltbRestraintPanel';
    div.innerHTML =
      '<div class="list" id="ltbRestraintList"></div>' +
      '<div class="addbar"><button type="button" id="addLtbRestraint">+ Add lateral restraint</button></div>' +
      '<div class="hint">The LTB boundary conditions of the two ends come from their degree-of-freedom flags under ' +
      'Geometry &amp; end conditions: U<sub>y</sub> = lateral displacement v held, R<sub>x</sub> = twist &phi; held (a fork end holds both), ' +
      'R<sub>z</sub> = lateral bending v&prime; held (laterally clamped end), warping = &phi;&prime; held. A free end holds none of them. ' +
      'Add intermediate restraints where purlins, ties or secondary beams hold the member. ' +
      'The LTB buckling length follows from these conditions and positions &mdash; the L<sub>E</sub> factor and the destabilising &times;1.2 switch ' +
      'do not affect EC3 LTB on the eigen route. Restraints that hold lateral displacement v also shorten the minor-axis strut length ' +
      'L<sub>cr,z</sub> to their spacing (SCI P360 6.2); the strut lengths otherwise follow the end fixities. ' +
      'A lateral cantilever (U<sub>y</sub> held at one end only) needs R<sub>z</sub> and R<sub>x</sub> at that end, or the lateral stiffness matrix is singular.</div>';
    host.parentNode.insertBefore(div, host.nextSibling);
    document.getElementById('addLtbRestraint').addEventListener('click', function () {
      if (!S.ltbRestraints) S.ltbRestraints = [];   // a state built without the field (an older fixture) - DEMO carries it since the 19 Sep 2026 review
      S.ltbRestraints.push({ pos: (S.L / 2).toFixed(3), v: true, phi: true, vp: false, phip: false });
      renderLtbRestraintList(); recompute();
    });
  }

  function renderLtbRestraintList() {
    var el = document.getElementById('ltbRestraintList');
    if (!el) return;
    el.innerHTML = (S.ltbRestraints || []).map(function (r, i) {
      return '<div class="row"><div class="rowhead">' +
        '<label style="flex:1">x, m <input type="number" step="0.01" data-i="' + i + '" data-k="pos" value="' + r.pos + '"></label>' +
        '<button type="button" class="del" data-del="' + i + '">Remove</button></div>' +
        '<div class="ltb-checks">' +
        chk(i, 'v', 'lateral v', r.v !== false) + chk(i, 'phi', 'twist &phi;', r.phi !== false) +
        chk(i, 'vp', 'lat. bending v&prime;', !!r.vp) + chk(i, 'phip', 'warping &phi;&prime;', !!r.phip) +
        '</div></div>';
    }).join('');
    function chk(i, k, lbl, on) {
      return '<label class="checkline"><input type="checkbox" data-i="' + i + '" data-k="' + k + '"' + (on ? ' checked' : '') + '> <span>' + lbl + '</span></label>';
    }
    el.querySelectorAll('input').forEach(function (inp) {
      inp.addEventListener('change', function () {
        var i = +inp.dataset.i, k = inp.dataset.k;
        S.ltbRestraints[i][k] = inp.type === 'checkbox' ? inp.checked : inp.value;
        recompute();
      });
    });
    el.querySelectorAll('button[data-del]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        S.ltbRestraints.splice(+btn.dataset.del, 1); renderLtbRestraintList(); recompute();
      });
    });
  }
  window.renderLtbRestraintList = renderLtbRestraintList;

  /* keep the panel in step with syncInputs(), and hide it when LTB is off */
  var _sync = window.syncInputs;
  window.syncInputs = function () {
    _sync.apply(this, arguments);
    injectUI(); renderLtbRestraintList();
    var show = (S.code === 'EC3' && (S.restraint || 'full') !== 'full');
    var p = document.getElementById('ltbRestraintPanel');
    if (p) p.style.display = show ? '' : 'none';
    var c1h = document.getElementById('c1Hint');
    if (c1h) c1h.innerHTML = (mcrMethod() === 'standard')
      ? 'Standard method: C<sub>1</sub> is derived from the moment diagram (NCCI SN003a tables for a simply supported UDL / central point load, ' +
        'the SCI end-moment curve for a linear gradient, otherwise the Serna quarter-point expression); it sets both M<sub>cr</sub> and ' +
        'k<sub>c</sub> = 1/&radic;C<sub>1</sub>. Override with a verified value (e.g. specialist software) if required.'
      : 'M<sub>cr</sub> is solved directly by the FE eigensolver &mdash; C<sub>1</sub> is not an input. ' +
        'It is back-calculated purely to form k<sub>c</sub> = 1/&radic;C<sub>1</sub> (NA 2.18). Override only to force k<sub>c</sub>; ' +
        'the eigen value is printed alongside.';
    /* the restraint panel is only meaningful to the eigen method; the standard
       route ignores intermediate restraints (whole-member LE = k x L) */
    var rp = document.getElementById('ltbRestraintPanel');
    if (rp && show) rp.style.opacity = (mcrMethod() === 'standard') ? '0.55' : '';
  };

  /* validate restraint positions */
  var _validate = window.validateInputs;
  window.validateInputs = function (py, E, uls, sls) {
    _validate.apply(this, arguments);
    (S.ltbRestraints || []).forEach(function (r, i) {
      var x = +r.pos;
      if (!isFinite(x) || x < 0 || x > S.L) throw 'Lateral restraint ' + (i + 1) + ' at x = ' + r.pos + ' m lies outside the span (0 to ' + S.L + ' m).';
    });
    if (S.Cmzo != null && !(isFinite(S.Cmzo) && S.Cmzo > 0)) throw 'C_mz override must be a positive number.';
  };

  assertSaggingPositive();
  if (typeof syncInputs === 'function') { try { injectUI(); renderLtbRestraintList(); } catch (e) {} }

  window.LTB_EIGEN = { mcrEigen: mcrEigen, solveLTB: solveLTB, momentFromSamples: momentFromSamples,
                       assertSaggingPositive: assertSaggingPositive };
})();
