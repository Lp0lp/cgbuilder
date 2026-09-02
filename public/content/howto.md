# How to use 🧱 CGBuilder

<div class="callout callout-note">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/></svg>
<p><strong>Before you start:</strong> bead-type prediction and the AA SMILES export both need a structure with <strong>explicit hydrogen atoms</strong> (e.g. a GROMACS <code>.gro</code> file, or a PDB run through a protonation tool). Everything else — loading, bead editing, .gro/.ndx/.map export — works regardless.</p>
</div>

<div class="toc">
<div class="toc-title">On this page</div>
<ol>
<li><a href="#step-load">Load a molecule</a></li>
<li><a href="#step-create">Create and select beads</a></li>
<li><a href="#step-assign">Assign atoms to a bead</a></li>
<li><a href="#step-edit">Edit bead properties</a></li>
<li><a href="#step-predict">Predict bead types</a></li>
<li><a href="#step-monitor">Monitor the mapping</a></li>
<li><a href="#step-sasa">SASA comparison</a></li>
<li><a href="#step-inspect">Inspect visually</a></li>
<li><a href="#step-export">Export your mapping</a></li>
</ol>
</div>

<h2 id="step-load"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg> 1. Load a molecule</h2>

Use the **Load molecule** panel in the top bar to bring in a structure:

- Click **Choose File** and select a `PDB`, or `GRO` file. The structure appears in the 3D viewer immediately.
- Click **Load example** to load a pre-built molecule with a complete Martini 3 mapping — useful for exploring the tool before working on your own system.
- Click **Load mapping** (enabled after loading a molecule) to paste an existing [Shaker](https://github.com/Lp0lp/shaker)-format mapping. The beads, types, and atom assignments will be applied automatically. **This replaces your current beads entirely** — unlike *Clear beads*, it does not ask for confirmation first.

<div class="callout callout-warning">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
<p>If the loaded structure has more than one residue, a warning appears above the bead list. This is usually a sign you've loaded more than you meant to (e.g. a molecule plus crystallographic waters, or a whole chain instead of a single fragment) — check it's really the molecule you intend to map.</p>
</div>

<h2 id="step-create"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><circle cx="12" cy="12" r="3"/></svg> 2. Create and select beads</h2>

- Click **New bead** to create a bead. The new bead becomes the active selection (highlighted in blue).
- Click a bead card in the list to select it. Click it again, or click empty space in the viewer, to deselect.
- Click **Clear beads** to remove all beads and start over (a confirmation dialog will appear).

<h2 id="step-assign"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 3l7.07 16.97 2.51-7.39 7.39-2.51L3 3z"/><path d="M13 13l6 6"/></svg> 3. Assign atoms to a bead</h2>

With a bead selected (highlighted in blue in the list):

- **Click an atom** in the 3D viewer to add it to the active bead. The bead sphere is placed at the geometric centre of its atoms.
- **Click the same atom again** to increase its weight, pulling the bead centre toward that atom. Useful when you want to tweak the placement of your CG bead.
- **Shift+click an atom** to reduce its weight by one step, or remove it entirely if the weight reaches zero.
- Atoms assigned to a bead are listed on its card and shown in the viewer.

<div class="callout callout-warning">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
<p>If a bead's boundary cuts through a bond at a heteroatom (N/O/S/P), a warning names the affected bead(s) above the bead list. See the Prediction module limitation below for more context. Consider keeping the whole functional group in one bead where possible.</p>
</div>

<h2 id="step-edit"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M11 4H4a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14a2 2 0 0 0 2-2v-7"/><path d="M18.5 2.5a2.121 2.121 0 0 1 3 3L12 15l-4 1 1-4 9.5-9.5z"/></svg> 4. Edit bead properties</h2>

Each bead card has three editable fields:

- **Name** — a unique label for this bead (e.g. `BB`, `SC1`).
- **Type** — the Martini 3 bead type (e.g. `P2`, `TC5`). Bead spheres and labels are colour-coded by polarity class:
  <span class="color-dot color-dot-c"></span>grey (C) &nbsp;
  <span class="color-dot color-dot-n"></span>blue (N) &nbsp;
  <span class="color-dot color-dot-p"></span>red (P) &nbsp;
  <span class="color-dot color-dot-q"></span>amber (Q/D) &nbsp;
  <span class="color-dot color-dot-x"></span>green (X)
- **Charge** — formal charge in units of *e*; most beads are 0.

<h2 id="step-predict"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="13 2 3 14 12 14 11 22 21 10 12 10 13 2"/></svg> 5. Predict bead types</h2>

- Click **Predict** on an individual bead card to get a type suggestion for that bead alone.
- Click **Predict bead types** (top of the bead list) to predict all beads at once.
- Predictions are inspired by the [AutoMartini](https://doi.org/10.1039/C9ME00183B) approach, using fragment free energies of transfer to suggest bead types. The implementation is independent and results may differ. A suggestion chip (`→ SP2a`) appears next to the Type field — click it to apply the suggestion.
- If the chemistry perceived from the structure implies a different formal charge than the bead's current **Charge** field (e.g. a deprotonated carboxylate detected), a matching suggestion chip appears next to **Charge** too. Applying it is always a separate, explicit click — predictions never silently change a value you haven't approved.

<div class="callout callout-warning">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
<p><strong>Always verify predictions manually.</strong> Automated assignments are a starting point — useful as initial guesses — not a final answer. Assignment algorithms are not perfect and can easily fail. Check each bead type against the chemical context of your molecule and the <a href="https://doi.org/10.1038/s41592-021-01098-3" target="_blank">Martini 3 guidelines</a> before using the mapping in a simulation.</p>
</div>

<div class="callout callout-note">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/></svg>
<p><strong>Limitation:</strong> predictions work best for chemically self-contained fragments. Ether oxygens that span two beads are capped with hydrogen during fragment construction, making them appear as hydroxyl groups — the predicted type will reflect hydroxyl character rather than ether. This is the sort of situation the capped-heteroatom warning above flags — manual correction of the bead type is recommended whenever it appears.</p>
</div>

If you use the bead type predictions, please also cite AutoMartini: <a href="https://doi.org/10.1021/acs.jctc.5c01178" target="_blank">doi:10.1021/acs.jctc.5c01178</a>.

<p class="info-footnote">Szczuka, Magdalena, et al. "Fast Parametrization of Martini3 Models for Fragments and Small Molecules." <em>Journal of Chemical Theory and Computation</em> 22.1 (2025): 610-623.</p>

<h2 id="step-monitor"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><line x1="18" y1="20" x2="18" y2="10"/><line x1="12" y1="20" x2="12" y2="4"/><line x1="6" y1="20" x2="6" y2="14"/></svg> 6. Monitor the mapping</h2>

The **Mapping** panel (top bar) tracks the quality of your mapping in real time. The mapping and mismatch calculations are only shown if bead types are defined for all beads. If any bead type is not set a warning is shown.

- **Heavy atoms** — total non-hydrogen atoms in the loaded structure. Atoms from period 4 or higher (Br, Se, I, ...) count as 2, since they're bulkier than a typical 2nd/3rd-period atom.
- **Beads** — number of beads defined, with a breakdown by size class (R = regular, S = small, T = tiny, U = virtual).
- **Mismatch** — difference between the heavy atoms your beads account for and the number expected from their size classes.
  <span class="color-dot color-dot-good"></span>Green = within tolerance &nbsp;
  <span class="color-dot color-dot-warn"></span>amber = acceptable &nbsp;
  <span class="color-dot color-dot-bad"></span>red = under- or over-mapped

<div class="callout callout-note">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/></svg>
<p>As per the Martini 3 guidelines, a mismatch of <strong>±1 per 10 heavy atoms</strong> is considered acceptable.</p>
</div>

<h2 id="step-sasa"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="12 2 2 7 12 12 22 7 12 2"/><polyline points="2 17 12 22 22 17"/><polyline points="2 12 12 17 22 12"/></svg> 7. SASA comparison</h2>

The **SASA** panel computes solvent-accessible surface areas for both representations, so you can check whether your mapping preserves the molecule's overall volume — this matters because volume directly affects interaction density in a simulation. 

The Å² <strong>numbers</strong> shown in this panel and the translucent <strong>surface mesh</strong> you see when toggling AA/CG surface come from two different, independent code paths. The numbers are computed by an implementation of the <a href="https://doi.org/10.1016/0022-2836(73)90011-9" target="_blank">Shrake-Rupley algorithm</a> (Shrake &amp; Rupley, 1973). The surface mesh is generated separately, by NGL's own surface algorithm, using the same 1.91&nbsp;Å probe radius for visual consistency but a different underlying method. 

- Toggle **AA surface** or **CG surface** to visualise the surfaces in the viewer.
- **Δ (CG−AA)** is the difference between the two Å² values.

<div class="callout callout-note">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/></svg>
<p>As per the Martini 3 guidelines, a Δ within <strong>±10%</strong> of the AA SASA is considered acceptable.</p>
</div>

<div class="callout callout-warning">
<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
<p>Absolute SASA values can differ slightly from other tools (e.g. GROMACS's <code>gmx sasa</code>) — SASA estimates are sensitive to the exact algorithm and no two implementations match perfectly. What matters is that AA and CG values are both computed with the <em>same</em> method so <strong>Δ stays a meaningful, internally-consistent comparison</strong>.</p>
</div>

<h2 id="step-inspect"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"/><circle cx="12" cy="12" r="3"/></svg> 8. Inspect visually</h2>

Use the **Viewer** panel toggles to interact with the molecular viewer:

- **AA labels** — show atom names on the all-atom structure.
- **CG labels** — show bead names on the CG spheres.
- **Solid beads** — make CG beads fully opaque (default is semi-transparent so the underlying atoms remain visible).
- **Light BG** — switch the viewport background between black and white. (Independent of the navbar's dark/light theme toggle, which only affects the page around the viewer.)
- **Save PNG** — export a 2× resolution screenshot of the current view.
- **Reset view** — re-centre and re-focus the camera.

<h2 id="step-export"><svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg> 9. Export your mapping</h2>

The output tabs at the bottom of the App page update live as you work. Use *Download* or *Copy*.

- **Shaker** — Python dict format used by [Shaker](https://github.com/Lp0lp/shaker).
- **Bartender** — plain-text mapping for [Bartender](https://github.com/Martini-Force-Field-Initiative/Bartender): one line per bead (`BEADS` section header, then `<bead number> <atom indices…>`, 1-based, atoms repeated by weight).
- **PyCGTOOL** — `.map` format for [PyCGTOOL](https://github.com/jag1g13/pycgtool): a `[ resname ]` section header followed by one line per bead (`name type charge atom1 atom2 …`), atoms repeated by weight.
- **.gro** — GROMACS coordinate file with one entry per bead at its mapped centre.
- **.ndx** — GROMACS index file grouping atoms by bead.
- **.map** — Martini .map file.
- **AA SMILES** — canonical SMILES string for the whole loaded structure, derived from its perceived chemistry. Needs explicit hydrogens (see the note at the top of this page); shows a guidance message instead if the structure has more than one disconnected fragment.

---

This version of CGBuilder is adapted from the original [CGBuilder by Jonathan Barnoud](https://github.com/jbarnoud/cgbuilder).
