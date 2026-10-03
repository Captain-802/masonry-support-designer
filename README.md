# Masonry Support Designer

Design check for a steel masonry support carried by a torsion beam:

- **BS 5950-1:2000** with SCI P057 (torsion), SCI P110 (biaxial stresses in the plate) and SCI P157 (stainless
  masonry support angles), or
- **Eurocode 3** with the UK National Annexes, with SCI P385 for the torsion beam.

## What it covers

- **Support member:** steel plate welded under the beam; hot-rolled angle, welded or bolted; cold-formed stainless
  angle; or the angle on its own, fixed to a wall or frame.
- **Loads:** masonry leaves (width, height and density), slab reactions, or direct dead and imposed line loads
  (any number of each).
- **Torsion beam:** UKPFC, UKB, UKC, RHS and SHS from the section list, or a user-defined channel, I-section or RHS.
- **Checks:** heel stress, biaxial stress ratio, toe deflection, welds and bolts; beam shear (including high shear),
  moment, lateral-torsional buckling, combined bending and torsion, local capacity and shear stress; twist and
  deflection at SLS. Optional construction-stage check for the plate.
- **Output:** pass/fail summary with utilisations, a cross-section drawn to scale, a full calculation sheet (Letter or
  A4) with formulas and substituted values, and engineering notes on the simplifications in the method.

## Validation

The BS 5950 route is benchmarked against 195 runs of an established commercial program for the same calculation:
worked examples, verification runs that drive every check to failure, and batches covering every support type,
channel, I and hollow beams, high shear, warping free and fixed, destabilising load, both twist limits and the
construction-stage check. The page recalculates all of them when it opens (Validation tab): 9,669 values, largest
difference 0.0002 %, and every pass/fail verdict agrees. The calculation sheet reproduces the benchmark program's
printed sheet line for line.

```
node tests/validate.js
```

The Eurocode 3 route has no external benchmark for the masonry support itself. Its beam checks come from a separately
benchmarked beam-design engine; the plate, angle, weld and bolt checks are verified by hand calculation.

## Files

| file | purpose |
|---|---|
| `index.html` | the page (built, do not edit by hand) |
| `src/template.html` | user interface: HTML, CSS and script |
| `src/sms.js` | BS 5950 calculation engine (N, mm) and calculation-sheet layout |
| `src/ec3_adapter.js` | Eurocode 3 route: loads and support checks around the beam engine |
| `src/beam-v03-engine.js`, `src/v03.css` | beam-design engine used by the Eurocode 3 route |
| `src/data.json` | section list and the 195 benchmark runs |
| `build.py` | inlines the engine and data into `index.html` |
| `tests/validate.js` | engine against the benchmark runs |

The page is one self-contained file: no external scripts, fonts or network calls, no `eval` and no inline event
handlers, so it can be embedded in a web page.

## Notes

Results are for checking by a competent engineer, who remains responsible for the design.
