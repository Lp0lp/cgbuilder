# 🧱 CGBuilder

**A browser-based visual editor for building Martini 3 coarse-grained molecule mappings.**

No installation is required. Open the site in a WebGL-capable browser:

**[Open CGBuilder](https://lp0lp.github.io/cgbuilder/)**

Load an all-atom structure, map beads interactively in a 3D viewer, get automated bead-type predictions, compare SASA values against your AA reference, and export ready-to-use mapping files.

> Fork of [jbarnoud/cgbuilder](https://github.com/jbarnoud/cgbuilder), substantially extended with weighted atom assignment, bead-type prediction, SASA comparison, [Shaker](https://github.com/Lp0lp/shaker)-format import/export, and more...

## Features

- **3D viewer** — interactive NGL.js viewport with licorice AA representation and coloured CG bead spheres overlaid in real time.
- **Weighted atom assignment** — click an atom once to add it to a bead; click again to increase its weight (pulling the bead centre toward that atom); Shift+click to reduce weight or remove.
- **Bead-type prediction** — inspired by [AutoMartini](https://doi.org/10.1039/C9ME00183B), using RDKit fragment free-energy-of-transfer heuristics. Suggestion chips appear next to each field; applying them is always an explicit click.
- **SASA comparison** — built-in Shrake-Rupley solver computes per-bead SASA and compares it against the all-atom reference surface rendered by NGL.
- **Live validation** — multi-residue warnings, capped-heteroatom bond-cutting warnings, atom-overlap detection, and bead-count/size guidelines displayed as you work.
- **[Shaker](https://github.com/Lp0lp/shaker)-format import** — paste an existing mapping to restore beads, types, charges, and atom assignments automatically.
- **Export** — `.gro` (CG coordinates), `.ndx` (atom-index groups), `.map` (martinize/backward mapping), [Shaker](https://github.com/Lp0lp/shaker) Python dict, [Bartender](https://github.com/Martini-Force-Field-Initiative/Bartender) mapping, and AA SMILES.

## Input formats

| Format | Load molecule | Notes                                                     |
| ------ | ------------- | --------------------------------------------------------- |
| `.pdb` | ✓             | Bond orders inferred; explicit H needed for prediction    |
| `.gro` | ✓             | GROMACS coordinate file; explicit H needed for prediction |

> **Explicit hydrogens required** for bead-type prediction and AA SMILES export. Use a protonation tool before loading if needed.

## Export formats

| File              | Contents                                                                                                                       |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `.gro`            | CG bead coordinates in GROMACS format                                                                                          |
| `.ndx`            | GROMACS index file — one group per bead, atoms by index                                                                        |
| `.map`            | Martini backward/martinize mapping — `[ to ]`/`[ martini ]`/`[ atoms ]` sections                                               |
| Shaker dict       | [Shaker](https://github.com/Lp0lp/shaker)-format Python assignment dict — bead names, types, charges, atom lists               |
| Bartender mapping | `BEADS` section + one line per bead with 1-based atom indices (repeated by weight)                                             |
| PyCGTOOL mapping  | [PyCGTOOL](https://github.com/jag1g13/pycgtool) `.map` — `[ resname ]` section + one line per bead (`name type charge atoms…`) |
| AA SMILES         | SMILES string for each bead's fragment (requires explicit H)                                                                   |

## Running tests

```bash
npm test             # single run
npm run test:watch   # watch mode
```

Tests cover all modules. No browser or DOM dependency is required; the suite runs in Node via Vitest with minimal stubs for NGL and RDKit globals.

---

## Dependencies

| Library                                           | Purpose                                                           |
| ------------------------------------------------- | ----------------------------------------------------------------- |
| [NGL](https://github.com/nglviewer/ngl)           | 3D molecular viewer (npm, bundled)                                |
| [@rdkit/rdkit](https://github.com/rdkit/rdkit-js) | SMILES generation and fragment canonicalisation (npm, code-split) |
| [marked](https://github.com/markedjs/marked)      | Markdown rendering for in-app docs (npm, bundled)                          |

Development tooling includes [TypeScript](https://www.typescriptlang.org/), [ESLint](https://eslint.org/), [Prettier](https://prettier.io/), [esbuild](https://esbuild.github.io/), [Vitest](https://vitest.dev/), [TypeDoc](https://typedoc.org/), [Husky](https://typicode.github.io/husky/), and [lint-staged](https://github.com/lint-staged/lint-staged).

## Self-deployment

CGBuilder is a static web application. You can build it locally and publish the generated `dist/` directory to any static hosting provider or web server.

Node.js 24 is required for building. The repository includes `.nvmrc` for version managers such as `nvm`. If you do not use `nvm`, install Node.js 24 through your preferred method.

```bash
git clone https://github.com/Lp0lp/cgbuilder.git
cd cgbuilder
nvm use
npm ci
npm run check
```

`npm run check` runs the lint, formatting, test, and production-build checks. After it completes, publish the contents of `dist/` to your hosting provider.

## Developer workflow

Start the development server with live reload:

```bash
npm run dev
```

Useful commands:

```bash
npm run fix           # apply ESLint fixes, then format files
npm run check         # lint, format check, tests, and build
npm test              # run tests once
npm run test:watch    # run tests in watch mode
npm run docs          # generate API documentation
```

- Code quality is enforced with [ESLint](https://eslint.org/).
- Formatting is enforced with [Prettier](https://prettier.io/).
- [Husky](https://typicode.github.io/husky/) automatically runs lint-staged on staged files before each commit. ESLint runs first, followed by Prettier.
- The GitHub Actions workflow repeats the checks in read-only mode for pull requests and validated builds.

## Acknowledgements

- Original CGBuilder by [@jbarnoud](https://github.com/jbarnoud/cgbuilder)
- Bead-type prediction approach inspired by [AutoMartini](https://doi.org/10.1021/acs.jctc.5b00056) (Bereau & Kremer, 2015) and extended for Martini 3 by [Szczuka et al., JCTC 2025](https://doi.org/10.1021/acs.jctc.5c01178)
- Shrake-Rupley SASA algorithm: [Shrake & Rupley, J. Mol. Biol. 1973](<https://doi.org/10.1016/0022-2836(73)90011-9>)
- [Martini 3 force field](http://cgmartini.nl/)
