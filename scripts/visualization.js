import { buildCanonTable, determineBeadType } from './prediction.js';
import {
    perceiveChemistry, fragmentToSmiles, moleculeToSmiles, beadDonorCount, countResidues,
    heavyAtomWeight, weightedHeavyAtomCount, cappedHeteroatoms,
} from './chemistry.js';
import { PROBE_RADIUS, aaSASA, cgSASA, beadsToPDB } from './sasa.js';
import { generateNDX, generateMap, generateGRO, generatePythonAssignments,
         download, copyTextToClipboard, bondAwareRepresentationParams,
         parseShakerMapping } from './fileformats.js';
import { loadRDKit } from './rdkit.js';

/* ===========================================================================
   Visualization — the main UI controller
   ===========================================================================
   Ties the bead-editing UI (the bead list panel, its fields and buttons) to
   the 3D viewer (NGL) and to the other modules: chemistry.js for bond order/
   charge perception, prediction.js for bead-type suggestions, sasa.js for
   the AA/CG surface area comparison, fileformats.js for import/export. Most
   of this file is DOM wiring and NGL representation management rather than
   novel algorithm — the one piece of real decision logic living here is
   _predictOneBead (the table-vs-Crippen lookup strategy for one bead),
   which belongs conceptually with prediction.js/chemistry.js even though it
   lives in this file. */

/**
 * RGB colour (0-1 range, for NGL) for a bead type's Martini polarity class,
 * read from the type code's first letter (or second, for a size-prefixed
 * type like "SP2a" / "TC5"). Matches the bead-card colour-coding described
 * in howto.md.
 * @param {string} beadType - Martini type code, e.g. "SP2a"
 * @returns {[number,number,number]}
 */
function typeColor(beadType) {
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
 * @param {Element} element
 * @param {string} className
 * @returns {Element|null}
 */
function findParentWithClass(element, className) {
    let node = element;
    while (node) {
        if (node.classList.contains(className)) return node;
        node = node.parentElement;
    }
    return null;
}

/**
 * Wire a tab's Copy button to copy the text of its output <pre>, with the
 * same "Copied!"/"Failed" button-text feedback used across all output tabs.
 * @param {string} buttonId
 * @param {string} outputId - id of the <pre> holding the text to copy
 */
function wireCopyButton(buttonId, outputId) {
    const button = document.getElementById(buttonId);
    if (!button) return;
    button.onclick = async () => {
        const originalText = button.textContent;
        const text = document.getElementById(outputId).textContent || "";
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
 * One instance is created per loaded molecule (see main.js).
 */
export class Visualization {
    /**
     * Wires up every static UI control that exists before a molecule is
     * loaded (label/CG/surface toggles, the Predict button, the output
     * tabs' Download/Copy buttons). Per-molecule state (the NGL component,
     * cached chemistry/SASA values) is populated later by attachRepresentation,
     * once a structure actually exists.
     * @param {object} collection - BeadCollection for this molecule
     * @param {object} stage - NGL Stage
     */
    constructor(collection, stage) {
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

        let toggleCGLabels = document.getElementById('toggle-cg-labels');
        toggleCGLabels.onchange = (e) => { this.showCGLabels = e.target.checked; this.drawCG(); };
        toggleCGLabels.checked = false;
        toggleCGLabels.disabled = false;

        let toggleCG = document.getElementById('toggle-cg');
        toggleCG.onchange = (event) => this.onToggleCG(event);
        toggleCG.checked = false;
        toggleCG.disabled = false;

        let toggleAASurface = document.getElementById('toggle-aa-surface');
        if (toggleAASurface) {
            toggleAASurface.onchange = (event) => this.onToggleAASurface(event);
            toggleAASurface.checked = false;
            toggleAASurface.disabled = false;
        }

        let toggleCGSurface = document.getElementById('toggle-cg-surface');
        if (toggleCGSurface) {
            toggleCGSurface.onchange = (event) => this.onToggleCGSurface(event);
            toggleCGSurface.checked = false;
            toggleCGSurface.disabled = false;
        }

        const predictBtn = document.getElementById('predict-types');
        if (predictBtn) {
            predictBtn.onclick = () => this.onPredictTypes();
            predictBtn.disabled = false;
        }

        document.getElementById('dl-ndx').onclick = () =>
            download('cgbuilder.ndx', generateNDX(this.collection));
        document.getElementById('dl-map').onclick = () =>
            download('cgbuilder.map', generateMap(this.collection));
        document.getElementById('dl-gro').onclick = () =>
            download('cgbuilder.gro', generateGRO(this.collection));
        document.getElementById('dl-py').onclick = () =>
            download('cgbuilder_assignments.py', generatePythonAssignments(this.collection));

        document.getElementById('dl-smiles').onclick = () =>
            download('cgbuilder.smi', document.getElementById('smiles-output').textContent || "");

        wireCopyButton('copy-py', 'py-output');
        wireCopyButton('copy-gro', 'gro-output');
        wireCopyButton('copy-ndx', 'ndx-output');
        wireCopyButton('copy-map', 'map-output');
        wireCopyButton('copy-smiles', 'smiles-output');
    }

    /** Shortcut for `this.collection.currentBead`. */
	get currentBead() {
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
     * @param {object} component - NGL StructureComponent for the loaded molecule
     */
    attachRepresentation(component) {
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
        // explicit hydrogens — see chemistry.js's module comment for why.
        // Computed once per load and cached, rather than per Predict click,
        // since it's also used to reflect real bond orders in the viewer.
        this.chemistry = perceiveChemistry(component.structure);

        const predictBtn = document.getElementById('predict-types');
        if (predictBtn) {
            predictBtn.disabled = !this.chemistry.available;
            predictBtn.title = this.chemistry.available
                ? ''
                : 'Bead-type prediction needs a structure with explicit hydrogen atoms '
                + '(e.g. a GROMACS .gro file, or a PDB run through a protonation tool).';
        }

        const residueWarning = document.getElementById('multi-residue-warning');
        if (residueWarning) residueWarning.hidden = countResidues(component.structure) <= 1;

        const smilesEl = document.getElementById('smiles-output');
        if (smilesEl) {
            if (!this.chemistry.available) {
                smilesEl.textContent = 'Needs a structure with explicit hydrogen atoms '
                    + '(e.g. a GROMACS .gro file, or a PDB run through a protonation tool).';
            } else {
                const smiles = moleculeToSmiles(component.structure, this.chemistry);
                smilesEl.textContent = smiles
                    ?? 'Could not generate a single SMILES — structure has more than one disconnected fragment.';
            }
        }

        if (this.chemistry.available) this._reflectBondOrders(component);
    }

    /**
     * Write perceiveChemistry's resolved bond orders back onto the actual
     * NGL structure (BondProxy.bondOrder is directly settable) and rebuild
     * representations so the viewer's ball+stick rendering shows real
     * double/triple bonds instead of whatever NGL guessed from a bare
     * PDB/GRO load (typically all-single).
     * @param {object} component - NGL StructureComponent
     */
    _reflectBondOrders(component) {
        const structure = component.structure;
        structure.eachAtom((atom) => {
            if (typeof atom.eachBond !== 'function') return;
            atom.eachBond((bond) => {
                const key = `${Math.min(bond.atomIndex1, bond.atomIndex2)}-${Math.max(bond.atomIndex1, bond.atomIndex2)}`;
                const order = this.chemistry.bondOrders.get(key);
                if (order !== undefined) bond.bondOrder = order;
            });
        });
        component.rebuildRepresentations();
    }

    /**
     * Show/hide the duplicate-atom-name warning for a just-loaded
     * structure. Duplicate names are a real problem since generateMap and
     * generatePythonAssignments (fileformats.js) both group/reference atoms
     * by name, not index.
     * @param {object} structure - NGL-style structure (eachAtom)
     */
    checkAtomNameUniqueness(structure) {
        const seen = new Set();
        let hasDupes = false;
        structure.eachAtom(ap => {
            const name = this.collection.atomName(ap);
            if (seen.has(name)) hasDupes = true;
            else seen.add(name);
        });
        document.getElementById('atom-name-warning').hidden = !hasDupes;
    }

    /**
     * Cache the structure's weighted heavy-atom count for updateMappingStats.
     * Weighted, not a plain count — period>=4 atoms (Br/Se/I...) count as 2,
     * so this stays on the same scale as bead-type-implied expected counts
     * (a bead correctly sized small around a bromine shouldn't look
     * "over-mapped" just because the raw atom count is low).
     * @param {object} structure - NGL-style structure (eachAtom)
     */
    countHeavyAtoms(structure) {
        this.nHeavyAtoms = weightedHeavyAtomCount(structure);
    }

    /**
     * Add a (initially hidden) per-atom label representation showing the
     * structure's original atom names (see BeadCollection.structureAtomNames),
     * and wire its visibility toggle.
     * @param {object} component - NGL StructureComponent
     */
    attachAALabels(component) {
        this.aa_labels = component.addRepresentation("label", {
            labelType: "text",
            labelText: this.collection.structureAtomNames(component.structure),
            labelGrouping: "atom",
            visible: false,
        });

        const toggle = document.getElementById('toggle-aa-labels');
        if (toggle) {
            toggle.checked = false;
            toggle.disabled = false;
            toggle.onchange = (event) => this.onToggleAALabels(event);
        }
    }

    /** "Solid beads" toggle: opacity only, see drawCG. */
    onToggleCG(event) {
        this.showCG = event.target.checked;
        this.drawCG();
    }

    /** "AA surface" toggle: just flips visibility, the representation already exists. */
    onToggleAASurface(event) {
        this.showAASurface = event.target.checked;
        if (this.aaSurface) this.aaSurface.setVisibility(this.showAASurface);
    }

    /** "CG surface" toggle: the representation is synthetic and rebuilt on demand, see drawCGSurface. */
    onToggleCGSurface(event) {
        this.showCGSurface = event.target.checked;
        this.drawCGSurface();
    }

    /**
     * "Predict bead types" button handler: load RDKit (cached after first
     * use, see rdkit.js), build the canonicalized lookup table once, then
     * run _predictOneBead for every bead in the collection.
     */
    async onPredictTypes() {
        const btn = document.getElementById('predict-types');
        if (!this.chemistry || !this.chemistry.available) return;
        const originalText = btn.textContent;
        btn.disabled = true;

        try {
            btn.textContent = 'Loading RDKit…';
            const RDKit = await loadRDKit();

            const canonTable = buildCanonTable(RDKit);

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
     * prediction.js's port of AutoMartini's algorithm meets a real bead):
     *
     *   1. Build the bead's fragment SMILES (chemistry.js's fragmentToSmiles,
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
     * @param {object} bead - a Bead with at least one atom assigned
     * @param {object} RDKit - the loaded RDKit_minimal module
     * @param {Object<string,number>} canonTable - buildCanonTable's result
     */
    _predictOneBead(bead, RDKit, canonTable) {
        if (bead.atoms.length === 0) return;
        const chemistry = this.chemistry;

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
        // chemistry.js's heavyAtomWeight. A ring or branch point (>=3 heavy
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
        let tableVal;
        if (inRing) {
            const keys = new Set();
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

        let deltaF, source;
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
     * flow as onPredictTypes, just for a single bead.
     * @param {object} bead
     * @param {Element} btn - the bead's own Predict button (for in-place feedback)
     */
    async onPredictBeadType(bead, btn) {
        if (!this.chemistry || !this.chemistry.available) return;
        const originalText = btn.textContent;
        btn.disabled = true;
        btn.textContent = '…';
        try {
            const RDKit = await loadRDKit();
            const canonTable = buildCanonTable(RDKit);
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
     * (sasa.js's beadsToPDB — the B-factor-as-radius trick, see its own
     * comment) and loads that into NGL as a new component carrying just
     * the surface representation. `_cgSurfaceToken` guards against a race
     * where the bead mapping changes again while the async `loadFile` is
     * still in flight — a stale load completing after a newer call has
     * already started (or the toggle's been switched off) is discarded
     * instead of replacing the current surface.
     */
    drawCGSurface() {
        if (this.cgSurfaceComp != null) {
            this.stage.removeComponent(this.cgSurfaceComp);
            this.cgSurfaceComp = null;
        }
        if (!this.showCGSurface) return;
        let pdb = beadsToPDB(this.collection);
        if (!pdb) return;

        let token = ++this._cgSurfaceToken;
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
    onToggleAALabels(event) {
        this.aa_labels.setVisibility(event.target.checked);
    }

    /**
     * NGL stage click handler (see stage.signals.clicked in main.js):
     * clicking an atom adds it to the current bead (or, with shift held,
     * removes it); clicking empty space deselects the current bead. A
     * click that hits an atom while no bead is selected is a no-op rather
     * than implicitly selecting one.
     * @param {object} pickingProxy - NGL PickingProxy for the click, or
     *   null/falsy for a click that hit nothing
     */
    onClick(pickingProxy) {
        if (pickingProxy && pickingProxy.atom) {
            if (!this.currentBead) return;
            if (pickingProxy.mouse && pickingProxy.mouse.shiftKey) {
                this.currentBead.removeAtom(pickingProxy.atom);
            } else {
                this.currentBead.addAtom(pickingProxy.atom);
            }
            this.updateSelection();
        } else if (!pickingProxy) {
            this.collection.deselectBead();
            this.updateSelection();
        }
    }

    /** "New bead" button handler — the new bead becomes selected (see BeadCollection.newBead). */
	onNewBead(event) {
	    this.collection.newBead();
	    this.updateSelection();
	}

    /**
     * Index of `node` within #bead-list's children — matches its position
     * in this.collection.beads, since createBeadList always rebuilds the
     * two in lockstep.
     * @param {Node} node
     * @returns {number} -1 if not found
     */
    _beadIndexForNode(node) {
        const nodes = document.getElementById("bead-list").childNodes;
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
     * @param {Event} event
     */
    onBeadSelected(event) {
        const tag = event.target.tagName;
        if (tag === "INPUT" || tag === "BUTTON" || tag === "FORM" || tag === "LABEL") return;

        const realTarget = findParentWithClass(event.target, "bead-view");
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
     * @param {Event} event
     */
	onBeadRemove(event) {
        const realTarget = findParentWithClass(event.target, "bead-view");
        const selected = this._beadIndexForNode(realTarget);
        if (selected >= 0) {
            this.collection.removeBead(selected);
            if (this.collection.beads.length === 0) this.collection.newBead();
            if (realTarget.classList.contains('selected-bead')) this.collection.selectBead(0);
        }
        this.updateSelection();
    }

    /**
     * Name field input handler: updates the bead's name as the user types,
     * then refreshes every output that includes bead names (updateName)
     * and re-checks for new name collisions (checkDuplicateNames).
     * @param {Event} event
     */
    onNameChange(event) {
        const realTarget = findParentWithClass(event.target, "bead-view");
        const index = this._beadIndexForNode(realTarget);
        if (index >= 0) this.collection.beads[index].name = event.target.value;
        this.updateName();
        this.checkDuplicateNames();
    }

    /**
     * NGL selection-language string for highlighting a bead's atoms in the
     * viewer: "@i,j,k" (atom indices, NGL's by-index selection syntax) for
     * a bead with atoms, or "not all" (selects nothing) otherwise.
     * @param {object} bead
     * @returns {string}
     */
	selectionString(bead) {
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
     */
    updateName() {
        this.updateNDX();
        this.updateMap();
        this.updateGRO();
        this.updatePY();
        this.updateSASA();
        this.updateMappingStats();
        this.updateCappedHeteroatomWarning();
        this.drawCG();
    }

    /**
     * Show/hide the "bead boundary cuts through a heteroatom" warning,
     * naming every affected bead. Beads whose boundary cuts through a
     * heteroatom (N/O/S/P) bond — capping that bond with hydrogen during
     * fragment construction can make a real ether/amine/thioether look
     * like a different, more terminal group than it actually is (see
     * chemistry.js's cappedHeteroatoms). Purely structural, so this is
     * shown regardless of chemistry.available.
     */
    updateCappedHeteroatomWarning() {
        const el = document.getElementById('capped-heteroatom-warning');
        if (!el) return;

        const names = [];
        for (const bead of this.collection.beads) {
            if (bead.atoms.length === 0) continue;
            if (cappedHeteroatoms(bead.atoms).length > 0) names.push(bead.name);
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
    updateMappingStats() {
        const heavyEl    = document.getElementById('map-heavy');
        const beadsEl    = document.getElementById('map-beads');
        const mismatchEl = document.getElementById('map-mismatch');
        if (!heavyEl) return;

        const reset = () => {
            [heavyEl, beadsEl, mismatchEl].forEach(el => { el.textContent = '—'; el.className = ''; });
        };

        if (!this.nHeavyAtoms) { reset(); return; }

        const nHeavy = this.nHeavyAtoms;
        const beads  = this.collection.beads;
        const nBeads = beads.length;
        heavyEl.textContent = nHeavy;

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
        const parts = [];
        if (counts.R) parts.push(`${counts.R}R`);
        if (counts.S) parts.push(`${counts.S}S`);
        if (counts.T) parts.push(`${counts.T}T`);
        if (counts.U) parts.push(`${counts.U}U`);
        beadsEl.textContent = `${nBeads} (${parts.join(' ')})`;

        const diff      = nHeavy - expectedHeavy;
        const tolerance = Math.max(1, Math.round(nHeavy / 10));
        const absDiff   = Math.abs(diff);
        const sign      = diff > 0 ? '+' : '';

        let label, cls;
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
    updateSelection() {
        let selString = this.selectionString(this.currentBead);
        this.representation.setSelection(selString);
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
     * @param {object} bead
     * @param {boolean} [isDuplicate] - whether this bead's name collides
     *   with another bead's (see createBeadList); shows the name field in
     *   its error state if so
     */
    createBeadListItem(bead, isDuplicate = false) {
        let list = document.getElementById("bead-list");
        let item = document.createElement("li");
        item.classList.add("bead-view");

        let headerRow = document.createElement("div");
        headerRow.classList.add("bead-header");

        let fieldsNode = document.createElement("div");
        fieldsNode.classList.add("bead-fields");

        const addLabeledField = (labelText, inputEl) => {
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
        let nameNode = document.createElement("input");
        nameNode.type = "text";
        nameNode.value = bead.name;
        nameNode.classList.add("bead-name");
        if (isDuplicate) nameNode.classList.add('input-error');
        nameNode.oninput = (event) => this.onNameChange(event);
        nameNode.addEventListener("mousedown", e => e.stopPropagation());
        addLabeledField("Name", nameNode);

        // TYPE
        let typeNode = document.createElement("input");
        typeNode.type = "text";
        typeNode.value = bead.type;
        typeNode.classList.add("bead-type");
        typeNode.oninput = (event) => { bead.type = event.target.value; this.updateName(); };
        typeNode.addEventListener("mousedown", e => e.stopPropagation());

        const typeWrap = document.createElement("div");
        typeWrap.classList.add("type-field-wrap");
        if (bead.suggestedType) {
            const chip = document.createElement("button");
            chip.classList.add("bead-type-chip");
            chip.textContent = `→ ${bead.suggestedType}`;
            chip.title = "Click to apply suggested type";
            chip.addEventListener("mousedown", e => e.stopPropagation());
            chip.onclick = (e) => {
                e.stopPropagation();
                bead.type = bead.suggestedType;
                typeNode.value = bead.suggestedType;
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
        let chargeNode = document.createElement("input");
        chargeNode.type = "number";
        chargeNode.step = "0.01";
        chargeNode.value = bead.charge;
        chargeNode.classList.add("bead-charge");
        chargeNode.oninput = (event) => { bead.charge = event.target.value; this.updateName(); };
        chargeNode.addEventListener("mousedown", e => e.stopPropagation());

        const chargeWrap = document.createElement("div");
        chargeWrap.classList.add("type-field-wrap");
        if (bead.suggestedCharge != null) {
            const chargeChip = document.createElement("button");
            chargeChip.classList.add("bead-type-chip");
            chargeChip.textContent = `→ ${bead.suggestedCharge}`;
            chargeChip.title = "Click to apply suggested charge";
            chargeChip.addEventListener("mousedown", e => e.stopPropagation());
            chargeChip.onclick = (e) => {
                e.stopPropagation();
                bead.charge = bead.suggestedCharge;
                chargeNode.value = bead.suggestedCharge;
                this.updateName();
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

        let removeNode = document.createElement("button");
        removeNode.textContent = "Delete";
        removeNode.classList.add("delete-bead", "btn-danger");
        removeNode.onclick = (event) => { event.stopPropagation(); this.onBeadRemove(event); };

        let predictNode = document.createElement("button");
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

            let nameList = document.createElement("ul");
            for (let i = 0; i < bead.atoms.length; i++) {
                const atom = bead.atoms[i];
                const name = this.collection.atomName(atom);
                const w = (bead.atomWeights && bead.atomWeights[atom.index])
                    ? bead.atomWeights[atom.index] : 1;
                let subitem = document.createElement("li");
                subitem.appendChild(document.createTextNode(w > 1 ? `${name} ×${w}` : name));
                if (this.collection.countBeadsForAtom(atom) > 1) {
                    let shareitem = document.createElement("abbr");
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
    createBeadList() {
        const counts = new Map();
        for (const bead of this.collection.beads)
            counts.set(bead.name, (counts.get(bead.name) || 0) + 1);
        const dupes = new Set([...counts.entries()].filter(([, n]) => n > 1).map(([k]) => k));
        for (const bead of this.collection.beads) this.createBeadListItem(bead, dupes.has(bead.name));
    }

    /**
     * Lighter-weight alternative to a full createBeadList rebuild: just
     * toggles each existing name field's error styling in place, for live
     * feedback while typing (called from onNameChange) without throwing
     * away and recreating every bead card on each keystroke.
     */
    checkDuplicateNames() {
        const counts = new Map();
        for (const bead of this.collection.beads)
            counts.set(bead.name, (counts.get(bead.name) || 0) + 1);
        const items = document.getElementById('bead-list').childNodes;
        let i = 0;
        for (const item of items) {
            const bead = this.collection.beads[i++];
            if (!bead) break;
            const nameInput = item.querySelector('.bead-name');
            if (nameInput) nameInput.classList.toggle('input-error', (counts.get(bead.name) || 0) > 1);
        }
    }

    /** Remove every bead card from the DOM (in preparation for createBeadList). */
    clearBeadList() {
        let list = document.getElementById('bead-list');
        while (list.lastChild) list.removeChild(list.lastChild);
    }

    /** Refresh the .ndx output tab. */
    updateNDX() { document.getElementById('ndx-output').textContent = generateNDX(this.collection); }
    /** Refresh the .map output tab. */
    updateMap() { document.getElementById('map-output').textContent = generateMap(this.collection); }
    /** Refresh the .gro output tab. */
    updateGRO() { document.getElementById('gro-output').textContent = generateGRO(this.collection); }
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
    updatePY() {
        const counts = new Map();
        for (const bead of this.collection.beads)
            counts.set(bead.name, (counts.get(bead.name) || 0) + 1);
        const hasDupes = [...counts.values()].some(n => n > 1);
        document.getElementById('py-warning').hidden = !hasDupes;
        document.getElementById('py-output').textContent = generatePythonAssignments(this.collection);
    }

    /**
     * Refresh the SASA panel: the cached AA value (computed once at load,
     * see attachRepresentation), the CG value (recomputed fresh every call,
     * since the bead mapping changes constantly), and their percentage
     * difference, colour-coded good/warn/bad at ±5%/±10%.
     */
    updateSASA() {
        const aaEl   = document.getElementById('aa-sasa');
        const cgEl   = document.getElementById('cg-sasa');
        const diffEl = document.getElementById('sasa-diff');
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
     * @param {string} text - pasted Shaker `mapping = {...}` text
     */
    loadShakerMapping(text) {
        if (!this.component) return;
        const beadDefs = parseShakerMapping(text);
        if (!beadDefs.length) { console.warn('No beads found in mapping file'); return; }

        // Collect name → index during iteration (avoids NGL proxy-reuse issues)
        const nameToIndex = new Map();
        this.component.structure.eachAtom(ap => {
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
                    bead.addAtom(this.component.structure.getAtomProxy(idx));
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
     * sasa.js), an optional name label per bead, and the selected bead
     * picked out in a fixed highlight colour rather than its own type
     * colour. The whole shape component is removed and recreated rather
     * than updated in place — simpler than diffing, and cheap enough at
     * the scale of a CG mapping's bead count. The "Solid beads" toggle is
     * just opacity (1 vs 0.4), not a different representation. Also
     * triggers a CG surface rebuild, since the surface depends on the same
     * bead positions.
     */
    drawCG() {
        let selectedColor = [0.25, 0.84, 0.96];
        let opacity = this.showCG ? 1 : 0.4;

        if (this.shapeComp != null) this.stage.removeComponent(this.shapeComp);

        let shape = new NGL.Shape("shape", {disablePicking: true});
        for (let bead of this.collection.beads) {
            const color = bead === this.currentBead ? selectedColor : typeColor(bead.type);
            if (bead.atoms.length > 0) {
                const center = bead.center;
                shape.addSphere(center, color, 1.12, bead.name);
                if (this.showCGLabels) shape.addText(
                    [center.x, center.y + 1.8, center.z], color, 2.5, bead.name
                );
            }
        }
        this.shapeComp = this.stage.addComponentFromObject(shape);
        this.shapeComp.addRepresentation("buffer", {opacity: opacity});

        this.drawCGSurface();
    }
}
