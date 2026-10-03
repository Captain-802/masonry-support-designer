// Recalculates every embedded benchmark run with the engine and reports the differences.
//   node tests/validate.js
'use strict';
const path = require('path');
const SMS = require(path.join(__dirname, '..', 'src', 'sms.js'));
const DATA = require(path.join(__dirname, '..', 'src', 'data.json'));

let n = 0, within = 0, worst = 0, at = '', cls = 0, clsBad = 0;
for (const c of DATA.cases) {
  const o = SMS.calc(c.p);
  for (const [k, ref] of Object.entries(c.refs)) {
    if (k === 'Class') { cls++; if (String(o.Class).toLowerCase() !== String(ref).toLowerCase()) clsBad++; continue; }
    const e = Math.abs(o[k] - ref) / Math.max(Math.abs(ref), 1e-9);
    n++;
    if (e <= 0.02) within++;
    if (!(e <= worst)) { worst = e; at = `${c.id} ${k}`; }
  }
}
console.log(`${DATA.cases.length} benchmark runs, ${n} values: ${within} within 2 %, largest difference ` +
  `${(100 * worst).toFixed(5)} % (${at}); section class ${cls - clsBad}/${cls} agree`);
process.exit(worst < 1e-5 && clsBad === 0 ? 0 : 1);
