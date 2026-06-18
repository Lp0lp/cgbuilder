# How to use 🧱 CGBuilder

## 1. Load a molecule

Use the **Load molecule** panel in the top bar to bring in a structure:

- Click **Choose File** and select a `PDB`, or `GRO` file. The structure appears in the 3D viewer immediately.
- Click **Load example** to load a pre-built molecule with a complete Martini 3 mapping — useful for exploring the tool before working on your own system.
- Click **Load mapping** (enabled after loading a molecule) to paste an existing [Shaker](https://github.com/Lp0lp/shaker)-format mapping. The beads, types, and atom assignments will be applied automatically.

## 2. Create and select beads

- Click **New bead** to create a bead. The new bead becomes the active selection (highlighted in blue).
- Click a bead card in the list to select it. Click it again, or click empty space in the viewer, to deselect.
- Click **Clear beads** to remove all beads and start over (a confirmation dialog will appear).

## 3. Assign atoms to a bead

With a bead selected (highlighted in blue in the list):

- **Click an atom** in the 3D viewer to add it to the active bead. The bead sphere is placed at the geometric centre of its atoms.
- **Click the same atom again** to increase its weight, pulling the bead centre toward that atom. Useful when you want to tweak the placement of your CG bead.
- **Shift+click an atom** to reduce its weight by one step, or remove it entirely if the weight reaches zero.
- Atoms assigned to a bead are listed on its card and shown in the viewer.

## 4. Edit bead properties

Each bead card has three editable fields:

- **Name** — a unique label for this bead (e.g. `BB`, `SC1`).
- **Type** — the Martini 3 bead type (e.g. `P2`, `TC5`). Bead spheres and labels are colour-coded by polarity class: grey (C), blue (N), red (P), amber (Q/D), green (X).
- **Charge** — formal charge in units of *e*; most beads are 0.

## 5. Predict bead types

- Click **Predict** on an individual bead card to get a type suggestion for that bead alone.
- Click **Predict bead types** (top of the bead list) to predict all beads at once.
- Predictions are inspired by the [AutoMartini](https://doi.org/10.1039/C9ME00183B) approach, using fragment free energies of transfer to suggest bead types. The implementation is independent and results may differ. A suggestion chip appears next to the Type field — click it to apply the suggestion.

> **Always verify predictions manually.** Automated assignments are a starting point, not a final answer. Check each bead type against the chemical context of your molecule and the [Martini 3 guidelines](https://doi.org/10.1038/s41592-021-01098-3) before using the mapping in a simulation.

**Limitation:** predictions work best for chemically self-contained fragments. Ether oxygens that span two beads are capped with hydrogen during fragment construction, making them appear as hydroxyl groups — the predicted type will reflect hydroxyl character rather than ether. Manual correction is recommended in such cases.

If you use the bead type predictions, please also cite AutoMartini (see [Resources](#resources)).

## 6. Monitor the mapping

The **Mapping** panel (top bar) tracks the quality of your mapping in real time:

- **Heavy atoms** — total non-hydrogen atoms in the loaded structure.
- **Beads** — number of beads defined, with a breakdown by size class (R = regular, S = small, T = tiny, U = virtual).
- **Mismatch** — difference between the heavy atoms your beads account for and the number expected from their size classes. A mismatch of ±1 per 10 heavy atoms is considered acceptable in Martini. Green = within tolerance, amber = acceptable, red = under- or over-mapped. Shown as a warning if any bead type is unset.

## 7. SASA comparison

The **SASA** panel computes solvent-accessible surface areas for both representations:

- Toggle **AA surface** or **CG surface** to visualise the surfaces in the viewer.
- **Δ (CG−AA)** is the difference in Å². A key goal of Martini mapping is to reproduce the molecular volume, as this directly affects interaction density in simulations. Aim for a Δ within ±10% of the AA SASA. 

## 8. Inspect visually

Use the **Viewer** panel toggles to interact with the molecular viewer:

- **AA labels** — show atom names on the all-atom structure.
- **CG labels** — show bead names on the CG spheres.
- **Solid beads** — make CG beads fully opaque (default is semi-transparent so the underlying atoms remain visible).
- **Light BG** — switch the viewport background between black and white.
- **Save PNG** — export a 2× resolution screenshot of the current view.
- **Reset view** — re-centre and re-focus the camera.

## 9. Export your mapping

The output tabs at the bottom of the App page update live as you work:

- **Shaker** — Python dict format used by [Shaker](https://github.com/Lp0lp/shaker). Use *Download* or *Copy*.
- **.gro** — GROMACS coordinate file with one entry per bead at its mapped centre.
- **.ndx** — GROMACS index file grouping atoms by bead.
- **.map** — Martini .map file.

---

Adapted from the original [CGBuilder by Jonathan Barnoud](https://github.com/jbarnoud/cgbuilder).
