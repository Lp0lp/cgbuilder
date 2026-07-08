# 🧱 CGBuilder

**A browser-based visual editor for building Martini 3 coarse-grained molecule mappings.**

Load an all-atom structure, map beads interactively in a 3D viewer, get automated bead-type predictions, compare SASA values against your AA reference, and export ready-to-use mapping files.

> Fork of [jbarnoud/cgbuilder](https://github.com/jbarnoud/cgbuilder), substantially extended with weighted atom assignment, bead-type prediction, SASA comparison, [Shaker](https://github.com/Lp0lp/shaker)-format import/export, and more...

---

## Features

- **3D viewer** — interactive NGL.js viewport with licorice AA representation and coloured CG bead spheres overlaid in real time.
- **Weighted atom assignment** — click an atom once to add it to a bead; click again to increase its weight (pulling the bead centre toward that atom); Shift+click to reduce weight or remove.
- **Bead-type prediction** — bead type prediction inspired by [AutoMartini](https://doi.org/10.1039/C9ME00183B), using RDKit fragment free-energy-of-transfer heuristics. Suggestion chips appear next to each field; applying them is always an explicit click.
- **SASA comparison** — built-in Shrake-Rupley solver computes per-bead SASA and compares it against the all-atom reference surface rendered by NGL.
- **Live validation** — multi-residue warnings, capped-heteroatom bond-cutting warnings, atom-overlap detection, and bead-count/size guidelines displayed as you work.
- **[Shaker](https://github.com/Lp0lp/shaker)-format import** — paste an existing mapping to restore beads, types, charges, and atom assignments automatically.
- **Export** — `.gro` (CG coordinates), `.ndx` (atom-index groups), `.map` (martinize/backward mapping), [Shaker](https://github.com/Lp0lp/shaker) Python dict, and AA SMILES (for fragment databases).
- **No build step** — plain ES modules served from a static file server.

---

## Quick start

```bash
# Any static file server works — file:// won't, because ES modules require HTTP.
npx serve .          # or: python3 -m http.server 8080
```

Open `http://localhost:3000` (or whichever port `serve` picks). Click **Load example** to try a pre-built Martini 3 mapping right away.

---

## Input formats

| Format | Load molecule | Notes |
|--------|--------------|-------|
| `.pdb` | ✓ | Bond orders inferred; explicit H needed for prediction |
| `.gro` | ✓ | GROMACS coordinate file; explicit H needed for prediction |
| `.sdf` / `.mol2` | ✓ | Bond-order metadata preserved; best for double-bond rendering |

> **Explicit hydrogens required** for bead-type prediction and AA SMILES export. Use a protonation tool before loading if needed.

---

## Export formats

| File | Contents |
|------|---------|
| `.gro` | CG bead coordinates in GROMACS format |
| `.ndx` | GROMACS index file — one group per bead, atoms by index |
| `.map` | Martini backward/martinize mapping — `[ to ]`/`[ martini ]`/`[ atoms ]` sections |
| Shaker dict | [Shaker](https://github.com/Lp0lp/shaker)-format Python assignment dict — bead names, types, charges, atom lists |
| Bartender mapping | `BEADS` section + one line per bead with 1-based atom indices (repeated by weight) |
| AA SMILES | SMILES string for each bead's fragment (requires explicit H) |

---

## Project layout

```
index.html              entry point
styles/main.css         all styling, CSS variables, dark-mode tokens
scripts/
  main.js               app bootstrap, file loading, event wiring
  visualization.js      NGL stage wrapper, bead rendering, UI interaction
  bead.js               Bead / BeadCollection
  chemistry.js          bond perception, valence assignment, RDKit SMILES
  prediction.js         bead-type prediction heuristics
  sasa.js               Shrake-Rupley SASA solver 
  fileformats.js        PDB/GRO/SDF/MOL2 parsing, .gro/.ndx/.map export
  rdkit.js              lazy RDKit-WASM loader (CDN, cached promise)
  example.js            bundled example molecule (PDB + mapping)
content/
  howto.md              in-app How to use documentation
  guidelines.md         Martini 3 bead-type parameterization reference material
tests/                  Vitest unit tests (87 tests across 6 files)
```

---

## Running tests

```bash
npm test          # single run
npm run test:watch   # watch mode
```

Tests cover `bead.js`, `chemistry.js`, `fileformats.js`, `prediction.js`, `sasa.js`, and `visualization.js`. No browser or DOM dependency — the test suite runs in Node via Vitest with minimal stubs for NGL and RDKit globals.

---

## Dependencies

All runtime dependencies are CDN-loaded; nothing to `npm install` for production use.

| Library | Version | Purpose |
|---------|---------|---------|
| [NGL.js](https://github.com/nglviewer/ngl) | v2.0.0-dev.39 | 3D molecular viewer |
| [RDKit.js](https://github.com/rdkit/rdkit-js) | latest CDN | SMILES generation and fragment canonicalisation |
| [marked](https://github.com/markedjs/marked) | v15 | Markdown rendering for in-app docs |
| [Inter](https://fonts.google.com/specimen/Inter) | — | UI typeface (Google Fonts) |

`devDependencies`: [Vitest](https://vitest.dev/) for the test suite.

---

## Acknowledgements

- Original CGBuilder by [@jbarnoud](https://github.com/jbarnoud/cgbuilder)
- Bead-type prediction approach inspired by [AutoMartini](https://doi.org/10.1021/acs.jctc.5b00056) (Bereau & Kremer, 2015) and extended for Martini 3 by [Szczuka et al., JCTC 2025](https://doi.org/10.1021/acs.jctc.5c01178)
- Shrake-Rupley SASA algorithm: [Shrake & Rupley, J. Mol. Biol. 1973](https://doi.org/10.1016/0022-2836(73)90011-9)
- [Martini 3 force field](http://cgmartini.nl/)
