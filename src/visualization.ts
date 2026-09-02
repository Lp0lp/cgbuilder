import { buildCanonTable, determineBeadType } from './prediction.js';
import {
    perceiveChemistry, fragmentToSmiles, moleculeToSmiles, beadDonorCount, countResidues,
    heavyAtomWeight, weightedHeavyAtomCount, cappedHeteroatoms,
} from './chemistry.js';
import { PROBE_RADIUS, aaSASA, cgSASA, beadsToPDB } from './sasa.js';
import { generateNDX, generateMap, generateGRO, generatePythonAssignments,
         generateBartender, generatePyCGTOOL, download, copyTextToClipboard,
         bondAwareRepresentationParams, parseShakerMapping } from './fileformats.js';
import { loadRDKit } from './rdkit.js';
import { byId } from './dom.js';
import { NGL } from './ngl.js';
import type { Bead, BeadCollection } from './bead.js';
import type {
    Chemistry, Component, NglRepresentation, PickingProxy, RDKitModule,
    Stage, StructureComponent,
} from './types.js';

/* ===========================================================================
   Visualization — the main UI controller
   ===========================================================================
   Ties the bead-editing UI (the bead list panel, its fields and buttons) to
   the 3D viewer (NGL) and to the other modules: chemistry.ts for bond order/
   charge perception, prediction.ts for bead-type suggestions, sasa.ts for
   the AA/CG surface area comparison, fileformats.ts for import/export. Most
   of this file is DOM wiring and NGL representation management rather than
   novel algorithm — the one piece of real decision logic living here is
   _predictOneBead (the table-vs-Crippen lookup strategy for one bead),
   which belongs conceptually with prediction.ts/chemistry.ts even though it
   lives in this file. */

/**
 * RGB colour (0-1 range, for NGL) for a bead type's Martini polarity class,
 * read from the type code's first letter (or second, for a size-prefixed
 * type like "SP2a" / "TC5"). Matches the bead-card colour-coding described
 * in howto.md.
 * @param beadType - Martini type code, e.g. "SP2a"
 */
export function typeColor(beadType: string): [number, number, number] {
    const t = (beadType || '').toUpperCase();
    const cls = (t[0] === 'S' || t[0] === 'T') ? t[1] : t[0];
    switch (cls) {
        case 'C': return [0.55, 0.55, 0.55]; // grey   – apolar
        case 'N': return [0.29, 0.56, 0.85]; // blue   – intermediate
        case 'P': return [0.91, 0.30, 0.24]; // red    – polar
        case 'Q': return [0.95, 0.61, 0.07]; // amber  – charged
        case 'D': return [0.95, 0.61, 0.07]; // amber  – divalent (charged)
        case 'X': return [0.18, 0.80, 0.44]; // green  – halogen
        default:  return [0.75, 0.75, 0.75]; // light grey – unknown/placeholder
    }
}

/**
 * Walk up from `element` to find the nearest ancestor (or itself) with
 * `className`, since DOM event targets are often a child of the element
 * the handler actually cares about (e.g. a click on a bead card's input
 * should still resolve to the `.bead-view` <li> for that bead).
 * @param element
 * @param className
 */
export function findParentWithClass(element: Element | null, className: string): Element | null {
    let node: Element | null = element;
    while (node) {
        if (node.classList.contains(className)) return node;
        node = node.parentElement;
    }
    return null;
}

/**
 * Wire a tab's Copy button to copy the text of its output <pre>, with the
 * same "Copied!"/"Failed" button-text feedback used across all output tabs.
 * @param buttonId
 * @param outputId - id of the <pre> holding the text to copy
 */
function wireCopyButton(buttonId: string, outputId: string): void {
    const button = byId(buttonId);
    if (!button) return;
    button.onclick = async () => {
        const originalText = button.textContent;
        const text = byId(outputId).textContent || "";
        try {
            await copyTextToClipboard(text);
            button.textContent = "Copied!";
            button.classList.add("copied");
            setTimeout(() => {
                button.textContent = originalText;
                button.classList.remove("copied");
            }, 1200);
        } catch (err) {
            button.textContent = "Failed";
            setTimeout(() => { button.textContent = originalText; }, 1200);
        }
    };
}

/**
 * Main UI controller: owns the NGL stage/representations, the bead-list DOM,
 * and the cached chemistry/SASA state for the currently loaded structure.
 * One instance is created per loaded molecule (see main.ts).
 */

/** Data-driven table for the simple output tabs (generator → DOM → download).
 *  Adding a new format = one row here; no other code needs to change. */
const OUTPUT_TABS = [
    { outputId: 'ndx-output',       dlId: 'dl-ndx',       copyId: 'copy-ndx',       filename: 'cgbuilder.ndx',           generate: generateNDX },
    { outputId: 'map-output',       dlId: 'dl-map',       copyId: 'copy-map',       filename: 'cgbuilder.map',           generate: generateMap },
    { outputId: 'gro-output',       dlId: 'dl-gro',       copyId: 'copy-gro',       filename: 'cgbuilder.gro',           generate: generateGRO },
    { outputId: 'bartender-output', dlId: 'dl-bartender', copyId: 'copy-bartender', filename: 'cgbuilder.bartender',     generate: generateBartender },
    { outputId: 'pycgtool-output',  dlId: 'dl-pycgtool',  copyId: 'copy-pycgtool',  filename: 'cgbuilder_pycgtool.map',  generate: generatePyCGTOOL },
] as const;

export class Visualization {
    collection: BeadCollection;
    representation: NglRepresentation | null;
    stage: Stage;
    shapeComp: Component | null;
    showCG: boolean;
    showCGLabels: boolean;

    // Solvent-accessible surface area (SASA) state.
    component: StructureComponent | null;   // the loaded AA structure component
    aaSurface: NglRepresentation | null;    // NGL surface representation on the AA structure
    showAASurface: boolean;
    cgSurfaceComp: Component | null;        // synthetic bead component carrying the CG surface
    showCGSurface: boolean;
    private _cgSurfaceToken: number;        // guards against stale async surface loads
    private _aaSASAValue: number | null;    // cached; recomputed only on molecule load
    nHeavyAtoms: number;
    private _canonTable: Record<string, number> | null; // cached RDKit-canonicalized lookup table
    chemistry: Chemistry | null;
    aa_labels: NglRepresentation | null;

    /**
     * Wires up every static UI control that exists before a molecule is
     * loaded (label/CG/surface toggles, the Predict button, the output
     * tabs' Download/Copy buttons). Per-molecule state (the NGL component,
     * cached chemistry/SASA values) is populated later by attachRepresentation,
     * once a structure actually exists.
     * @param collection - BeadCollection for this molecule
     * @param stage - NGL Stage
     */
    constructor(collection: BeadCollection, stage: Stage) {
        this.collection = collection;
        this.representation = null;
        this.stage = stage;
        this.shapeComp = null;
        this.showCG = false;
        this.showCGLabels = false;

        // Solvent-accessible surface area (SASA) state.
        this.component = null;        // the loaded AA structure component
        this.aaSurface = null;        // NGL surface representation on the AA structure
        this.showAASurface = false;
        this.cgSurfaceComp = null;    // synthetic bead component carrying the CG surface
        this.showCGSurface = false;
        this._cgSurfaceToken = 0;     // guards against stale async surface loads
        this._aaSASAValue = null;     // cached; recomputed only on molecule load
        this.nHeavyAtoms = 0;
        this._canonTable = null;      // cached RDKit-canonicalized lookup table, see _loadPredictionDeps
        this.chemistry = null;
        this.aa_labels = null;

        const toggleCGLabels = byId<HTMLInputElement>('toggle-cg-labels');
        toggleCGLabels.onchange = (e) => { this.showCGLabels = (e.target as HTMLInputElement).checked; this.drawCG(); };
        toggleCGLabels.checked = false;
        toggleCGLabels.disabled = false;

        const toggleCG = byId<HTMLInputElement>('toggle-cg');
        toggleCG.onchange = (event) => this.onToggleCG(event);
        toggleCG.checked = false;
        toggleCG.disabled = false;

        const toggleAASurface = byId<HTMLInputElement>('toggle-aa-surface');
        if (toggleAASurface) {
            toggleAASurface.onchange = (event) => this.onToggleAASurface(event);
            toggleAASurface.checked = false;
            toggleAASurface.disabled = false;
        }

        const toggleCGSurface = byId<HTMLInputElement>('toggle-cg-surface');
        if (toggleCGSurface) {
            toggleCGSurface.onchange = (event) => this.onToggleCGSurface(event);
            toggleCGSurface.checked = false;
            toggleCGSurface.disabled = false;
        }

        const predictBtn = byId<HTMLButtonElement>('predict-types');
        if (predictBtn) {
            predictBtn.onclick = () => this.onPredictTypes();
            predictBtn.disabled = false;
        }

        for (const tab of OUTPUT_TABS) {
            const { dlId, copyId, outputId, filename, generate } = tab;
            byId<HTMLButtonElement>(dlId).onclick = () => download(filename, generate(this.collection));
            wireCopyButton(copyId, outputId);
        }
        byId<HTMLButtonElement>('dl-py').onclick = () =>
            download('cgbuilder_assignments.py', generatePythonAssignments(this.collection));
        byId<HTMLButtonElement>('dl-smiles').onclick = () =>
            download('cgbuilder.smi', byId('smiles-output').textContent || "");
        wireCopyButton('copy-py', 'py-output');
        wireCopyButton('copy-smiles', 'smiles-output');
    }

    /** Shortcut for `this.collection.currentBead`. */
    get currentBead(): Bead | null {
        return this.collection.currentBead;
    }

    /**
     * Set up everything that depends on a just-loaded AA structure: the
     * highlighted-selection ball+stick representation, the AA surface
     * representation (hidden until toggled), the cached whole-structure
     * SASA value, chemistry perception (gating the Predict button and the
     * AA SMILES tab on whether explicit hydrogens are present), the
     * multi-residue warning, and — only when chemistry is available —
     * writing real bond orders back onto the structure (_reflectBondOrders)
     * so the viewer shows real double/triple bonds instead of NGL's
     * geometry-based guess.
     * @param component - NGL StructureComponent for the loaded molecule
     */
    attachRepresentation(component: StructureComponent): void {
        this.component = component;
        this.representation = component.addRepresentation(
            "ball+stick",
            bondAwareRepresentationParams({
                sele: "not all",
                radiusScale: 1.6,
                color: "#f4b642",
                opacity: 0.6
            }),
        );

        this.aaSurface = component.addRepresentation("surface", {
            surfaceType: "sas",
            probeRadius: PROBE_RADIUS,
            color: "#f4b642",
            opacity: 0.6,
            wireframe: true,
            visible: this.showAASurface,
            useWorker: false,
        });

        this._aaSASAValue = aaSASA(component.structure, PROBE_RADIUS);
        this.updateSASA();

        // Chemistry perception (bond order / aromaticity / charge) requires
        // explicit hydrogens — see chemistry.ts's module comment for why.
        // Computed once per load and cached, rather than per Predict click,
        // since it's also used to reflect real bond orders in the viewer.
        const chem = perceiveChemistry(component.structure);
        this.chemistry = chem;

        const predictBtn = byId<HTMLButtonElement>('predict-types');
        if (predictBtn) {
            predictBtn.disabled = !chem.available;
            predictBtn.title = chem.available
                ? ''
                : 'Bead-type prediction needs a structure with explicit hydrogen atoms '
                + '(e.g. a GROMACS .gro file, or a PDB run through a protonation tool).';
        }

        const residueWarning = byId('multi-residue-warning');
        if (residueWarning) residueWarning.hidden = countResidues(component.structure) <= 1;

        const smilesEl = byId('smiles-output');
        if (smilesEl) {
            if (!chem.available) {
                smilesEl.textContent = 'Needs a structure with explicit hydrogen atoms '
                    + '(e.g. a GROMACS .gro file, or a PDB run through a protonation tool).';
            } else {
                const smiles = moleculeToSmiles(component.structure, chem);
                smilesEl.textContent = smiles
                    ?? 'Could not generate a single SMILES — structure has more than one disconnected fragment.';
            }
        }

        if (chem.available) this._reflectBondOrders(component);
    }

    /**
     * Write perceiveChemistry's resolved bond orders back onto the actual
     * NGL structure (BondProxy.bondOrder is directly settable) and rebuild
     * representations so the viewer's ball+stick rendering shows real
     * double/triple bonds instead of whatever NGL guessed from a bare
     * PDB/GRO load (typically all-single).
     * @param component - NGL StructureComponent
     */
    _reflectBondOrders(component: StructureComponent): void {
        const chem = this.chemistry;
        if (!chem) return;
        const structure = component.structure;
        structure.eachAtom((atom) => {
            if (typeof atom.eachBond !== 'function') return;
            atom.eachBond((bond) => {
                const key = `${Math.min(bond.atomIndex1, bond.atomIndex2)}-${Math.max(bond.atomIndex1, bond.atomIndex2)}`;
                const order = chem.bondOrders.get(key);
                if (order !== undefined) bond.bondOrder = order;
            });
        });
        component.rebuildRepresentations();
    }

    /**
     * Show/hide the duplicate-atom-name warning for a just-loaded
     * structure. Duplicate names are a real problem since generateMap and
     * generatePythonAssignments (fileformats.ts) both group/reference atoms
     * by name, not index.
     * @param structure - NGL-style structure (eachAtom)
     */
    checkAtomNameUniqueness(structure: StructureComponent['structure']): void {
        const seen = new Set<string>();
        let hasDupes = false;
        structure.eachAtom(ap => {
            const name = this.collection.atomName(ap);
            if (seen.has(name)) hasDupes = true;
            else seen.add(name);
        });
        byId('atom-name-warning').hidden = !hasDupes;
    }

    /**
     * Cache the structure's weighted heavy-atom count for updateMappingStats.
     * Weighted, not a plain count — period>=4 atoms (Br/Se/I...) count as 2,
     * so this stays on the same scale as bead-type-implied expected counts
     * (a bead correctly sized small around a bromine shouldn't look
     * "over-mapped" just because the raw atom count is low).
     * @param structure - NGL-style structure (eachAtom)
     */
    countHeavyAtoms(structure: StructureComponent['structure']): void {
        this.nHeavyAtoms = weightedHeavyAtomCount(structure);
    }

    /**
     * Add a (initially hidden) per-atom label representation showing the
     * structure's original atom names (see BeadCollection.structureAtomNames),
     * and wire its visibility toggle.
     * @param component - NGL StructureComponent
     */
    attachAALabels(component: StructureComponent): void {
        this.aa_labels = component.addRepresentation("label", {
            labelType: "text",
            labelText: this.collection.structureAtomNames(component.structure),
            labelGrouping: "atom",
            visible: false,
        });

        const toggle = byId<HTMLInputElement>('toggle-aa-labels');
        if (toggle) {
            toggle.checked = false;
            toggle.disabled = false;
            toggle.onchange = (event) => this.onToggleAALabels(event);
        }
    }

    /** "Solid beads" toggle: opacity only, see drawCG. */
    onToggleCG(event: Event): void {
        this.showCG = (event.target as HTMLInputElement).checked;
        this.drawCG();
    }

    /** "AA surface" toggle: just flips visibility, the representation already exists. */
    onToggleAASurface(event: Event): void {
        this.showAASurface = (event.target as HTMLInputElement).checked;
        if (this.aaSurface) this.aaSurface.setVisibility(this.showAASurface);
    }

    /** "CG surface" toggle: the representation is synthetic and rebuilt on demand, see drawCGSurface. */
    onToggleCGSurface(event: Event): void {
        this.showCGSurface = (event.target as HTMLInputElement).checked;
        this.drawCGSurface();
    }

    /**
     * Load RDKit (cached across the whole page session, see rdkit.ts) and
     * this molecule's RDKit-canonicalized lookup table (cached on `this`,
     * since buildCanonTable re-canonicalizes ~150 table entries through
     * RDKit's WASM module every time it's called — too expensive to redo on
     * every individual bead's "Predict" click).
     */
    async _loadPredictionDeps(): Promise<{ RDKit: RDKitModule; canonTable: Record<string, number> }> {
        const RDKit = await loadRDKit();
        if (!this._canonTable) this._canonTable = buildCanonTable(RDKit);
        return { RDKit, canonTable: this._canonTable };
    }

    /**
     * "Predict bead types" button handler: load RDKit and the canonicalized
     * lookup table (cached, see _loadPredictionDeps), then run
     * _predictOneBead for every bead in the collection.
     */
    async onPredictTypes(): Promise<void> {
        const btn = byId<HTMLButtonElement>('predict-types');
        if (!this.chemistry || !this.chemistry.available) return;
        const originalText = btn.textContent;
        btn.disabled = true;

        try {
            btn.textContent = 'Loading RDKit…';
            const { RDKit, canonTable } = await this._loadPredictionDeps();

            btn.textContent = 'Predicting…';
            for (const bead of this.collection.beads) {
                this._predictOneBead(bead, RDKit, canonTable);
            }

            this.updateSelection();
        } catch (err) {
            console.error('Bead type prediction failed:', err);
            btn.textContent = 'Prediction failed';
            setTimeout(() => {
                btn.textContent = originalText;
                btn.disabled = false;
            }, 2000);
            return;
        }

        btn.textContent = originalText;
        btn.disabled = false;
    }

    /**
     * Predict one bead's Martini type (and, separately, suggest a charge),
     * writing the results onto `bead.suggestedType`/`bead.suggestedCharge`
     * rather than returning them — the bead-card UI reads these directly to
     * show the suggestion chips. The actual decision logic (this is where
     * prediction.ts's port of AutoMartini's algorithm meets a real bead):
     *
     *   1. Build the bead's fragment SMILES (chemistry.ts's fragmentToSmiles,
     *      H-capped at the bead boundary) and get RDKit's canonical form
     *      plus descriptors (Crippen logP, H-bond acceptor count) from it.
     *   2. Derive charge/halogen-presence/ring-membership/weighted heavy
     *      count/ring-or-branch directly from the resolved chemistry —
     *      charge is never read from the bead's own (user-editable) charge
     *      field, only ever derived from the structure.
     *   3. Look up a reference deltaF: for a ring-containing bead, try
     *      every aromatic-notation serialization of the fragment (one per
     *      possible DFS start atom, since the table only stores one
     *      direction per fragment) against canonTable first; fall back to
     *      the plain canonical-SMILES table lookup; if nothing matches
     *      either way, estimate deltaF from the fragment's own Crippen
     *      logP instead (logged as the "crippen" source vs. "table").
     *   4. Call determineBeadType with all of the above to get the actual
     *      predicted type code.
     *
     * Also separately flags a suggested charge update (bead.suggestedCharge)
     * whenever the chemistry-derived charge differs from what's currently
     * in the bead's own charge field — e.g. a deprotonated carboxylate
     * detected from explicit hydrogens that the user hasn't reflected yet.
     * Only ever a suggestion: applying it is a separate, explicit UI action.
     * @param bead - a Bead with at least one atom assigned
     * @param RDKit - the loaded RDKit_minimal module
     * @param canonTable - buildCanonTable's result
     */
    _predictOneBead(bead: Bead, RDKit: RDKitModule, canonTable: Record<string, number>): void {
        if (bead.atoms.length === 0) return;
        const chemistry = this.chemistry;
        if (!chemistry) return;

        const smiles = fragmentToSmiles(bead.atoms, chemistry);
        if (!smiles) return;

        const mol = RDKit.get_mol(smiles);
        if (!mol) {
            console.warn(`Bead "${bead.name}": invalid SMILES "${smiles}" — skipping`);
            return;
        }

        const canonSmiles = mol.get_smiles();
        const desc = JSON.parse(mol.get_descriptors());
        mol.delete();

        const charge = bead.atoms.reduce((s, a) => s + (chemistry.charges.get(a.index) ?? 0), 0);
        const hasHalogen = bead.atoms.some(a =>
            ['F','CL','BR','I'].includes((a.element || '').toUpperCase()));
        const inRing     = bead.atoms.some(a => chemistry.aromaticAtoms.has(a.index));
        const heavyAtoms = bead.atoms.filter(a => (a.element || 'C').toUpperCase() !== 'H');
        // Period->=4 atoms (Br/Se/I...) count as 2 toward bead size — see
        // chemistry.ts's heavyAtomWeight. A ring or branch point (>=3 heavy
        // neighbours anywhere in the bead) downgrades a weighted count of 4
        // from R to S, per the Martini 3 SI's bead-size convention.
        const weightedHeavyCount = heavyAtoms.reduce((s, a) => s + heavyAtomWeight(a.element), 0);
        const ringOrBranched = bead.atoms.some(a =>
            chemistry.ringAtoms.has(a.index) || chemistry.branchAtoms.has(a.index));
        const hDonors    = beadDonorCount(bead, chemistry);

        // Table-first lookup. The AutoMartini table uses open-chain
        // aromatic SMILES ("cc", "cn", "ncs"…) for partial-ring fragments
        // — invalid for RDKit but kept verbatim in canonTable. For aromatic
        // fragments, the serialisation depends on the DFS start atom (a
        // C–N pair is "cn" from C but "nc" from N), and the table only
        // stores one direction, so we enumerate the aromatic-notation
        // SMILES from every start atom and try each. The RDKit canonical
        // (Kekulé) form is the final fallback so full rings (c1ccccc1) and
        // non-aromatic entries still match. Real double bonds (inRing=false)
        // skip the aromatic lookup and fall through to Crippen → TC4.
        let lookupKey = canonSmiles;
        let tableVal: number | undefined;
        if (inRing) {
            const keys = new Set<string>();
            for (const a of bead.atoms) {
                if (!chemistry.aromaticAtoms.has(a.index)) continue;
                const k = fragmentToSmiles(bead.atoms, chemistry, {
                    aromaticNotation: true, startIndex: a.index });
                if (k) keys.add(k);
            }
            for (const k of keys) {
                if (canonTable[k] !== undefined) { lookupKey = k; tableVal = canonTable[k]; break; }
            }
            if (tableVal === undefined && keys.size) lookupKey = [...keys][0];
        }
        if (tableVal === undefined) tableVal = canonTable[canonSmiles];

        let deltaF: number, source: string;
        if (tableVal !== undefined) {
            deltaF = tableVal;
            source = 'table';
        } else {
            const logP = desc.CrippenClogP ?? 0;
            deltaF = 0.008314 * 300 * Math.LN10 * logP;
            source = 'crippen';
        }

        const hAcceptors = desc.NumHBA ?? 0;
        console.log(
            `Bead "${bead.name}": ${smiles} → key=${lookupKey} | `
            + `δf=${deltaF.toFixed(1)} (${source}) `
            + `arom=${inRing} HBD=${hDonors} HBA=${hAcceptors}`);

        bead.suggestedType = determineBeadType(
            { deltaF, charge, hDonors, hAcceptors, hasHalogen, inRing, weightedHeavyCount, ringOrBranched });

        // Suggest updating the Charge field too when the bead's atoms imply
        // a charge the user hasn't already set — e.g. a deprotonated
        // carboxylate detected from explicit hydrogens. Only surfaced when
        // it would actually change something.
        bead.suggestedCharge = (charge !== 0 && charge !== (bead.charge || 0)) ? charge : null;
    }

    /**
     * Per-bead "Predict" button handler — same RDKit-load-then-predict
     * flow as onPredictTypes (see _loadPredictionDeps), just for a single
     * bead.
     * @param bead
     * @param btn - the bead's own Predict button (for in-place feedback)
     */
    async onPredictBeadType(bead: Bead, btn: HTMLButtonElement): Promise<void> {
        if (!this.chemistry || !this.chemistry.available) return;
        const originalText = btn.textContent;
        btn.disabled = true;
        btn.textContent = '…';
        try {
            const { RDKit, canonTable } = await this._loadPredictionDeps();
            this._predictOneBead(bead, RDKit, canonTable);
            this.updateSelection();
        } catch (err) {
            console.error('Per-bead prediction failed:', err);
            btn.textContent = 'Failed';
            setTimeout(() => { btn.textContent = originalText; btn.disabled = false; }, 2000);
            return;
        }
        btn.textContent = originalText;
        btn.disabled = false;
    }

    /**
     * (Re)build the CG surface representation: removes any existing one,
     * and, if the toggle is on, serializes the beads to a synthetic PDB
     * (sasa.ts's beadsToPDB — the B-factor-as-radius trick, see its own
     * comment) and loads that into NGL as a new component carrying just
     * the surface representation. `_cgSurfaceToken` guards against a race
     * where the bead mapping changes again while the async `loadFile` is
     * still in flight — a stale load completing after a newer call has
     * already started (or the toggle's been switched off) is discarded
     * instead of replacing the current surface.
     */
    drawCGSurface(): void {
        if (this.cgSurfaceComp != null) {
            this.stage.removeComponent(this.cgSurfaceComp);
            this.cgSurfaceComp = null;
        }
        if (!this.showCGSurface) return;
        const pdb = beadsToPDB(this.collection);
        if (!pdb) return;

        const token = ++this._cgSurfaceToken;
        this.stage
            .loadFile(new Blob([pdb], {type: "text/plain"}), {ext: "pdb"})
            .then((comp) => {
                if (token !== this._cgSurfaceToken || !this.showCGSurface) {
                    this.stage.removeComponent(comp);
                    return;
                }
                comp.addRepresentation("surface", {
                    surfaceType: "sas",
                    radiusType: "bfactor",
                    radiusScale: 1.0,
                    probeRadius: PROBE_RADIUS,
                    color: "#7fc8a9",
                    opacity: 0.6,
                    wireframe: true,
                    useWorker: false,
                });
                this.cgSurfaceComp = comp;
            })
            .catch((err) => console.error("Error building CG surface:", err));
    }

    /** "AA labels" toggle. */
    onToggleAALabels(event: Event): void {
        this.aa_labels?.setVisibility((event.target as HTMLInputElement).checked);
    }

    /**
     * NGL stage click handler (see stage.signals.clicked in main.ts):
     * clicking an atom adds it to the current bead (or, with shift held,
     * removes it); clicking empty space deselects the current bead. A
     * click that hits an atom while no bead is selected is a no-op rather
     * than implicitly selecting one.
     * @param pickingProxy - NGL PickingProxy for the click, or
     *   null/undefined for a click that hit nothing
     */
    onClick(pickingProxy: PickingProxy | undefined): void {
        if (pickingProxy && pickingProxy.atom) {
            const bead = this.currentBead;
            if (!bead) return;
            if (pickingProxy.mouse && pickingProxy.mouse.shiftKey) {
                bead.removeAtom(pickingProxy.atom);
            } else {
                bead.addAtom(pickingProxy.atom);
            }
            this.updateSelection();
        } else if (!pickingProxy) {
            this.collection.deselectBead();
            this.updateSelection();
        }
    }

    /** "New bead" button handler — the new bead becomes selected (see BeadCollection.newBead). */
    onNewBead(event: Event): void {
        this.collection.newBead();
        this.updateSelection();
    }

    /**
     * Index of `node` within #bead-list's children — matches its position
     * in this.collection.beads, since createBeadList always rebuilds the
     * two in lockstep.
     * @param node
     * @returns -1 if not found
     */
    _beadIndexForNode(node: Node | null): number {
        const nodes = byId("bead-list").childNodes;
        let index = 0;
        for (const child of nodes) {
            if (child === node) return index;
            index += 1;
        }
        return -1;
    }

    /**
     * Bead-card click handler: selects that bead, or deselects it if it was
     * already the current one (a toggle). Ignores clicks that land on an
     * actual form control (input/button/etc.) within the card, so editing
     * a field doesn't also change the selection.
     * @param event
     */
    onBeadSelected(event: Event): void {
        const tag = (event.target as HTMLElement).tagName;
        if (tag === "INPUT" || tag === "BUTTON" || tag === "FORM" || tag === "LABEL") return;

        const realTarget = findParentWithClass(event.target as Element, "bead-view");
        const index = this._beadIndexForNode(realTarget);
        if (index >= 0) {
            if (this.collection.beads[index] === this.currentBead) {
                this.collection.deselectBead();
            } else {
                this.collection.selectBead(index);
            }
        }
        this.updateSelection();
    }

    /**
     * "Delete" button handler for one bead card. Keeps the invariant that
     * there's always at least one bead (creates a fresh one if the
     * collection would otherwise be empty), and re-selects bead 0 if the
     * bead being removed was the selected one.
     * @param event
     */
    onBeadRemove(event: Event): void {
        const realTarget = findParentWithClass(event.target as Element, "bead-view");
        const selected = this._beadIndexForNode(realTarget);
        if (selected >= 0) {
            this.collection.removeBead(selected);
            if (this.collection.beads.length === 0) this.collection.newBead();
            if (realTarget!.classList.contains('selected-bead')) this.collection.selectBead(0);
        }
        this.updateSelection();
    }

    /**
     * Name field input handler: updates the bead's name as the user types,
     * then refreshes every output that includes bead names (updateName)
     * and re-checks for new name collisions (checkDuplicateNames).
     * @param event
     */
    onNameChange(event: Event): void {
        const realTarget = findParentWithClass(event.target as Element, "bead-view");
        const index = this._beadIndexForNode(realTarget);
        if (index >= 0) this.collection.beads[index].name = (event.target as HTMLInputElement).value;
        this.updateName(false); // a name can't affect bead position or R/S/T size class
        this.checkDuplicateNames();
    }

    /**
     * NGL selection-language string for highlighting a bead's atoms in the
     * viewer: "@i,j,k" (atom indices, NGL's by-index selection syntax) for
     * a bead with atoms, or "not all" (selects nothing) otherwise.
     * @param bead
     */
    selectionString(bead: Bead | null): string {
        if (bead && bead.atoms.length > 0) {
            let sel = "@";
            for (let i = 0; i < bead.atoms.length; i++) {
                if (sel !== '@') sel = sel + ',';
                sel = sel + bead.atoms[i].index;
            }
            return sel;
        }
        return "not all";
    }

    /**
     * Recompute every output that depends on the current bead mapping
     * (despite the name — historically just the bead name, now everything
     * downstream of it): the four export tabs, the SASA comparison, the
     * Mapping panel stats, the chopped-heteroatom warning, and the CG
     * viewer. Called after essentially any bead edit (name/type/charge
     * change, atom add/remove, bead add/remove/select).
     * @param rebuildSurface - forwarded to drawCG; pass false for edits
     *   (name, charge) that can't affect bead position or R/S/T size class,
     *   so the CG surface mesh isn't needlessly torn down and rebuilt. The
     *   export tabs/SASA/mapping stats above are unaffected either way,
     *   since they're cheap and always reflect the latest values.
     */
    updateName(rebuildSurface = true): void {
        this.updateOutputTabs();
        this.updatePY();
        this.updateSASA();
        this.updateMappingStats();
        this.updateCappedHeteroatomWarning();
        this.drawCG(rebuildSurface);
    }

    /**
     * Show/hide the "bead boundary cuts through a heteroatom" warning,
     * naming every affected bead. Beads whose boundary cuts through a
     * heteroatom (N/O/S/P) bond — capping that bond with hydrogen during
     * fragment construction can make a real ether/amine/thioether look
     * like a different, more terminal group than it actually is (see
     * chemistry.ts's cappedHeteroatoms). Purely structural, so this is
     * shown regardless of chemistry.available.
     */
    updateCappedHeteroatomWarning(): void {
        const el = byId('capped-heteroatom-warning');
        if (!el) return;

        const names: string[] = [];
        for (const bead of this.collection.beads) {
            if (bead.atoms.length === 0) continue;
            if (cappedHeteroatoms(bead.atoms).length > 0) names.push(bead.name ?? "");
        }

        if (names.length === 0) {
            el.hidden = true;
            return;
        }
        el.hidden = false;
        const plural = names.length > 1;
        el.textContent = `⚠ Bead${plural ? 's' : ''} ${names.join(', ')} `
            + `${plural ? 'have' : 'has'} a bond crossing the bead boundary at a heteroatom `
            + `(N/O/S/P) — predictions may reflect a more terminal group (e.g. an alcohol or `
            + `primary amine) than the real, larger functional group. Consider keeping the `
            + `whole group in one bead.`;
    }

    /**
     * Refresh the Mapping panel: heavy-atom count, bead count (with an
     * R/S/T/U size-class breakdown), and the mismatch indicator comparing
     * the structure's real weighted heavy-atom count against what the
     * beads' own sizes imply (R=4, S=3, T=2, U=0 heavy atoms each).
     * Requires every bead to have a real type assigned first (not the
     * "TYPe" placeholder) — shows an "Assign bead types" prompt instead of
     * a mismatch verdict until then. Mismatch tolerance is ±1 per 10 heavy
     * atoms (rounded, minimum 1): within tolerance is "Acceptable", exact
     * is "OK", otherwise under/over-mapped depending on sign.
     */
    updateMappingStats(): void {
        const heavyEl    = byId('map-heavy');
        const beadsEl    = byId('map-beads');
        const mismatchEl = byId('map-mismatch');
        if (!heavyEl) return;

        const reset = (): void => {
            [heavyEl, beadsEl, mismatchEl].forEach(el => { el.textContent = '—'; el.className = ''; });
        };

        if (!this.nHeavyAtoms) { reset(); return; }

        const nHeavy = this.nHeavyAtoms;
        const beads  = this.collection.beads;
        const nBeads = beads.length;
        heavyEl.textContent = String(nHeavy);

        const allTyped = nBeads > 0 && beads.every(b => b.type && b.type !== 'TYPe');

        if (!allTyped) {
            beadsEl.textContent = `${nBeads}`;
            mismatchEl.textContent = '⚠ Assign bead types';
            mismatchEl.className = 'sasa-diff-warn';
            return;
        }

        const counts = { R: 0, S: 0, T: 0, U: 0 };
        let expectedHeavy = 0;
        for (const bead of beads) {
            const t = (bead.type || '').toUpperCase();
            if      (t[0] === 'T') { counts.T++; expectedHeavy += 2; }
            else if (t[0] === 'S') { counts.S++; expectedHeavy += 3; }
            else if (t[0] === 'U') { counts.U++;                     }
            else                   { counts.R++; expectedHeavy += 4; }
        }
        const parts: string[] = [];
        if (counts.R) parts.push(`${counts.R}R`);
        if (counts.S) parts.push(`${counts.S}S`);
        if (counts.T) parts.push(`${counts.T}T`);
        if (counts.U) parts.push(`${counts.U}U`);
        beadsEl.textContent = `${nBeads} (${parts.join(' ')})`;

        const diff      = nHeavy - expectedHeavy;
        const tolerance = Math.max(1, Math.round(nHeavy / 10));
        const absDiff   = Math.abs(diff);
        const sign      = diff > 0 ? '+' : '';

        let label: string, cls: string;
        if      (absDiff === 0)        { label = 'OK';           cls = 'sasa-diff-good'; }
        else if (absDiff <= tolerance) { label = 'Acceptable';   cls = 'sasa-diff-warn'; }
        else if (diff > 0)             { label = 'Under-mapped'; cls = 'sasa-diff-bad';  }
        else                           { label = 'Over-mapped';  cls = 'sasa-diff-bad';  }

        mismatchEl.textContent = `${sign}${diff}  ${label}`;
        mismatchEl.className = cls;
    }

    /**
     * Full redraw after any change to which beads exist or which is
     * selected: re-highlights the current bead's atoms in the AA viewer,
     * rebuilds the entire bead-list DOM from scratch, and cascades into
     * updateName for every downstream output.
     */
    updateSelection(): void {
        const selString = this.selectionString(this.currentBead);
        this.representation?.setSelection(selString);
        this.clearBeadList();
        this.createBeadList();
        this.updateName();
    }

    /**
     * Build and append the DOM for one bead card: Name/Type/Charge fields
     * (each with a suggestion "chip" button when a prediction is available),
     * Delete/Predict buttons, and — if the bead has any atoms — a
     * collapsible atom list with a 🔗 indicator on any atom shared with
     * another bead. Applies the `selected-bead` highlight and scrolls it
     * into view if this is the current bead.
     * @param bead
     * @param isDuplicate - whether this bead's name collides with another
     *   bead's (see createBeadList); shows the name field in its error
     *   state if so
     */
    createBeadListItem(bead: Bead, isDuplicate = false): void {
        const list = byId("bead-list");
        const item = document.createElement("li");
        item.classList.add("bead-view");

        const headerRow = document.createElement("div");
        headerRow.classList.add("bead-header");

        const fieldsNode = document.createElement("div");
        fieldsNode.classList.add("bead-fields");

        const addLabeledField = (labelText: string, inputEl: HTMLElement): void => {
            const wrap = document.createElement("div");
            wrap.classList.add("field");
            const lab = document.createElement("div");
            lab.classList.add("field-label");
            lab.textContent = labelText;
            wrap.appendChild(lab);
            wrap.appendChild(inputEl);
            fieldsNode.appendChild(wrap);
        };

        // NAME
        const nameNode = document.createElement("input");
        nameNode.type = "text";
        nameNode.value = bead.name ?? "";
        nameNode.classList.add("bead-name");
        if (isDuplicate) nameNode.classList.add('input-error');
        nameNode.oninput = (event) => this.onNameChange(event);
        nameNode.addEventListener("mousedown", e => e.stopPropagation());
        addLabeledField("Name", nameNode);

        // TYPE
        const typeNode = document.createElement("input");
        typeNode.type = "text";
        typeNode.value = bead.type;
        typeNode.classList.add("bead-type");
        typeNode.oninput = (event) => { bead.type = (event.target as HTMLInputElement).value; this.updateName(); };
        typeNode.addEventListener("mousedown", e => e.stopPropagation());

        const typeWrap = document.createElement("div");
        typeWrap.classList.add("type-field-wrap");
        if (bead.suggestedType) {
            const suggestedType = bead.suggestedType;
            const chip = document.createElement("button");
            chip.classList.add("bead-type-chip");
            chip.textContent = `→ ${suggestedType}`;
            chip.title = "Click to apply suggested type";
            chip.addEventListener("mousedown", e => e.stopPropagation());
            chip.onclick = (e) => {
                e.stopPropagation();
                bead.type = suggestedType;
                typeNode.value = suggestedType;
                this.updateName();
            };
            typeWrap.appendChild(chip);
        }

        const typeFieldEl = document.createElement("div");
        typeFieldEl.classList.add("field");
        const typeLab = document.createElement("div");
        typeLab.classList.add("field-label");
        typeLab.textContent = "Type";
        typeFieldEl.appendChild(typeLab);
        typeFieldEl.appendChild(typeNode);
        typeFieldEl.appendChild(typeWrap);
        fieldsNode.appendChild(typeFieldEl);

        // CHARGE
        const chargeNode = document.createElement("input");
        chargeNode.type = "number";
        chargeNode.step = "0.01";
        chargeNode.value = String(bead.charge);
        chargeNode.classList.add("bead-charge");
        chargeNode.oninput = (event) => { bead.charge = (event.target as HTMLInputElement).value; this.updateName(false); };
        chargeNode.addEventListener("mousedown", e => e.stopPropagation());

        const chargeWrap = document.createElement("div");
        chargeWrap.classList.add("type-field-wrap");
        if (bead.suggestedCharge != null) {
            const suggestedCharge = bead.suggestedCharge;
            const chargeChip = document.createElement("button");
            chargeChip.classList.add("bead-type-chip");
            chargeChip.textContent = `→ ${suggestedCharge}`;
            chargeChip.title = "Click to apply suggested charge";
            chargeChip.addEventListener("mousedown", e => e.stopPropagation());
            chargeChip.onclick = (e) => {
                e.stopPropagation();
                bead.charge = suggestedCharge;
                chargeNode.value = String(suggestedCharge);
                this.updateName(false);
            };
            chargeWrap.appendChild(chargeChip);
        }

        const chargeFieldEl = document.createElement("div");
        chargeFieldEl.classList.add("field");
        const chargeLab = document.createElement("div");
        chargeLab.classList.add("field-label");
        chargeLab.textContent = "Charge";
        chargeFieldEl.appendChild(chargeLab);
        chargeFieldEl.appendChild(chargeNode);
        chargeFieldEl.appendChild(chargeWrap);
        fieldsNode.appendChild(chargeFieldEl);

        // DELETE + PREDICT button column
        const btnCol = document.createElement("div");
        btnCol.classList.add("bead-btn-col");

        const removeNode = document.createElement("button");
        removeNode.textContent = "Delete";
        removeNode.classList.add("delete-bead", "btn-danger");
        removeNode.onclick = (event) => { event.stopPropagation(); this.onBeadRemove(event); };

        const predictNode = document.createElement("button");
        predictNode.textContent = "Predict";
        predictNode.classList.add("predict-bead");
        predictNode.addEventListener("mousedown", e => e.stopPropagation());
        predictNode.onclick = (e) => { e.stopPropagation(); this.onPredictBeadType(bead, predictNode); };

        btnCol.appendChild(removeNode);
        btnCol.appendChild(predictNode);

        headerRow.appendChild(fieldsNode);
        headerRow.appendChild(btnCol);
        item.appendChild(headerRow);

        // ATOM LIST — collapsible
        if (bead.atoms.length > 0) {
            const atomDetails = document.createElement("details");
            atomDetails.classList.add("atom-list-details");

            const atomSummary = document.createElement("summary");
            atomSummary.classList.add("atom-list-summary");
            const n = bead.atoms.length;
            atomSummary.textContent = `${n} atom${n !== 1 ? 's' : ''}`;
            atomSummary.addEventListener("click", e => e.stopPropagation());
            atomDetails.appendChild(atomSummary);

            const nameList = document.createElement("ul");
            for (let i = 0; i < bead.atoms.length; i++) {
                const atom = bead.atoms[i];
                const name = this.collection.atomName(atom);
                const w = (bead.atomWeights && bead.atomWeights[atom.index])
                    ? bead.atomWeights[atom.index] : 1;
                const subitem = document.createElement("li");
                subitem.appendChild(document.createTextNode(w > 1 ? `${name} ×${w}` : name));
                if (this.collection.countBeadsForAtom(atom) > 1) {
                    const shareitem = document.createElement("abbr");
                    shareitem.title = "This atom is shared between multiple beads.";
                    shareitem.textContent = " 🔗";
                    subitem.appendChild(shareitem);
                }
                nameList.appendChild(subitem);
            }
            atomDetails.appendChild(nameList);
            item.appendChild(atomDetails);
        }

        item.onclick = (event) => this.onBeadSelected(event);
        list.appendChild(item);

        if (bead === this.currentBead) {
            item.classList.add("selected-bead");
            item.scrollIntoView({ block: 'nearest' });
        }
    }

    /**
     * Rebuild the whole bead-list DOM from the current collection,
     * detecting name collisions once up front and flagging every bead
     * involved (not just the second-or-later occurrence) via
     * createBeadListItem's isDuplicate flag.
     */
    createBeadList(): void {
        const counts = new Map<string, number>();
        for (const bead of this.collection.beads) {
            const name = bead.name ?? "";
            counts.set(name, (counts.get(name) || 0) + 1);
        }
        const dupes = new Set([...counts.entries()].filter(([, n]) => n > 1).map(([k]) => k));
        for (const bead of this.collection.beads) this.createBeadListItem(bead, dupes.has(bead.name ?? ""));
    }

    /**
     * Lighter-weight alternative to a full createBeadList rebuild: just
     * toggles each existing name field's error styling in place, for live
     * feedback while typing (called from onNameChange) without throwing
     * away and recreating every bead card on each keystroke.
     */
    checkDuplicateNames(): void {
        const counts = new Map<string, number>();
        for (const bead of this.collection.beads) {
            const name = bead.name ?? "";
            counts.set(name, (counts.get(name) || 0) + 1);
        }
        const items = byId('bead-list').childNodes;
        let i = 0;
        for (const item of items) {
            const bead = this.collection.beads[i++];
            if (!bead) break;
            const nameInput = (item as Element).querySelector<HTMLElement>('.bead-name');
            if (nameInput) nameInput.classList.toggle('input-error', (counts.get(bead.name ?? "") || 0) > 1);
        }
    }

    /** Remove every bead card from the DOM (in preparation for createBeadList). */
    clearBeadList(): void {
        const list = byId('bead-list');
        while (list.lastChild) list.removeChild(list.lastChild);
    }

    updateOutputTabs(): void {
        for (const tab of OUTPUT_TABS) {
            byId(tab.outputId).textContent = tab.generate(this.collection);
        }
    }
    /**
     * Refresh the Shaker output tab, including its own duplicate-bead-name
     * warning — generatePythonAssignments writes one dict-literal line per
     * bead, so duplicate names still both appear in the generated text, but
     * Python itself silently keeps only the last one when it evaluates a
     * dict literal with a repeated key. So the breakage happens downstream,
     * when Shaker reads the file, not in our own generation step — but
     * it's flagged here as an error (not a cosmetic duplicate) since the
     * output is genuinely wrong, just not wrong yet at the point we show it.
     */
    updatePY(): void {
        const counts = new Map<string, number>();
        for (const bead of this.collection.beads) {
            const name = bead.name ?? "";
            counts.set(name, (counts.get(name) || 0) + 1);
        }
        const hasDupes = [...counts.values()].some(n => n > 1);
        byId('py-warning').hidden = !hasDupes;
        byId('py-output').textContent = generatePythonAssignments(this.collection);
    }

    /**
     * Refresh the SASA panel: the cached AA value (computed once at load,
     * see attachRepresentation), the CG value (recomputed fresh every call,
     * since the bead mapping changes constantly), and their percentage
     * difference, colour-coded good/warn/bad at ±5%/±10%.
     */
    updateSASA(): void {
        const aaEl   = byId('aa-sasa');
        const cgEl   = byId('cg-sasa');
        const diffEl = byId('sasa-diff');
        if (!aaEl || !cgEl || !diffEl) return;

        const aaVal = this._aaSASAValue;
        aaEl.textContent = aaVal !== null ? aaVal.toFixed(1) : '—';

        const cgVal = cgSASA(this.collection, PROBE_RADIUS);
        cgEl.textContent = cgVal > 0 ? cgVal.toFixed(1) : '—';

        if (aaVal !== null && aaVal > 0 && cgVal > 0) {
            const pct = (cgVal - aaVal) / aaVal * 100;
            const sign = pct >= 0 ? '+' : '';
            diffEl.textContent = `${sign}${pct.toFixed(1)}%`;
            diffEl.className = Math.abs(pct) < 5  ? 'sasa-diff-good'
                             : Math.abs(pct) < 10 ? 'sasa-diff-warn'
                             :                      'sasa-diff-bad';
        } else {
            diffEl.textContent = '—';
            diffEl.className = '';
        }
    }

    /**
     * "Load mapping" feature: parse a pasted Shaker assignments dict
     * (parseShakerMapping) and replace the current bead collection with
     * it. Atoms are matched by name (via BeadCollection.atomName, so the
     * original source-file names if recorded), not index, since a pasted
     * mapping has no idea what index a name corresponds to in THIS
     * particular load of the structure. Any mapping atom name not found in
     * the structure is skipped with a console warning rather than failing
     * the whole import. No-ops if no molecule is loaded yet, or if the
     * text contains no recognizable bead definitions.
     * @param text - pasted Shaker `mapping = {...}` text
     */
    loadShakerMapping(text: string): void {
        const component = this.component;
        if (!component) return;
        const beadDefs = parseShakerMapping(text);
        if (!beadDefs.length) { console.warn('No beads found in mapping file'); return; }

        // Collect name → index during iteration (avoids NGL proxy-reuse issues)
        const nameToIndex = new Map<string, number>();
        component.structure.eachAtom(ap => {
            const name = this.collection.atomName(ap);
            if (!nameToIndex.has(name)) nameToIndex.set(name, ap.index);
        });

        this.collection.clearBeads();

        for (const def of beadDefs) {
            const bead = this.collection.newBead();
            bead.name = def.name;
            bead.type = def.type;
            bead.charge = def.charge;
            for (const atomName of def.atoms) {
                const idx = nameToIndex.get(atomName);
                if (idx !== undefined) {
                    bead.addAtom(component.structure.getAtomProxy(idx));
                } else {
                    console.warn(`Mapping import: atom "${atomName}" not found in structure`);
                }
            }
        }

        this.collection.selectBead(0);
        this.updateSelection();
    }

    /**
     * Redraw the CG representation from scratch: a fresh NGL.Shape with
     * one sphere per non-empty bead (radius is fixed/cosmetic here, unlike
     * the real per-size-class radius used for SASA/the CG surface — see
     * sasa.ts), an optional name label per bead, and the selected bead
     * picked out in a fixed highlight colour rather than its own type
     * colour. The whole shape component is removed and recreated rather
     * than updated in place — simpler than diffing, and cheap enough at
     * the scale of a CG mapping's bead count. The "Solid beads" toggle is
     * just opacity (1 vs 0.4), not a different representation. Also
     * triggers a CG surface rebuild by default, since the surface depends
     * on the same bead positions.
     * @param rebuildSurface - pass false to skip the trailing drawCGSurface
     *   call, for edits that can't affect bead position or R/S/T size class
     *   (see updateName) — the surface mesh has no other dependency, so
     *   rebuilding it would be wasted, async, main-thread work.
     */
    drawCG(rebuildSurface = true): void {
        const selectedColor: [number, number, number] = [0.25, 0.84, 0.96];
        const opacity = this.showCG ? 1 : 0.4;

        if (this.shapeComp != null) this.stage.removeComponent(this.shapeComp);

        const shape = new NGL.Shape("shape", {disablePicking: true});
        for (const bead of this.collection.beads) {
            const color = bead === this.currentBead ? selectedColor : typeColor(bead.type);
            if (bead.atoms.length > 0) {
                const center = bead.center;
                shape.addSphere(center, color, 1.12, bead.name ?? "");
                if (this.showCGLabels) shape.addText(
                    [center.x, center.y + 1.8, center.z], color, 2.5, bead.name ?? ""
                );
            }
        }
        this.shapeComp = this.stage.addComponentFromObject(shape);
        this.shapeComp.addRepresentation("buffer", {opacity: opacity});

        if (rebuildSurface) this.drawCGSurface();
    }
}
