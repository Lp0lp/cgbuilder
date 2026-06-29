/* ===========================================================================
   Bead / BeadCollection — the CG mapping model
   ===========================================================================
   A Bead is one coarse-grained particle: a name, a Martini bead type, a
   formal charge, and a weighted set of all-atom atoms mapped onto it. A
   BeadCollection is the full mapping for one structure — every bead, which
   one is currently selected for editing, and a record of the structure's
   original atom names (see setOriginalAtomNames) so exported files can use
   the real source-file names rather than whatever NGL parsed them as. */

/**
 * One coarse-grained bead: name, Martini type, formal charge, and a
 * weighted set of all-atom atoms.
 */
export class Bead {
	constructor () {
		this._name = null;
        this._type = "TYPe";
		this._charge = 0;
		this.atoms = [];
        this.atomWeights = {};        // key: atom.index -> integer weight
	}

    /**
     * Position of `atom` in this bead's atoms array.
     * @param {object} atom - atom proxy (has an `index`)
     * @returns {number} array index, or -1 if not in this bead
     */
    indexOf(atom) {
        for (let i = 0; i < this.atoms.length; i++) {
            if (this.atoms[i].index === atom.index) return i;
        }
        return -1;
    }

    /**
     * Add an atom to this bead, or increment its weight by 1 if it's
     * already in the bead — repeated clicks on the same atom pull the
     * bead's centre (see `center`) further toward it without duplicating
     * the atom in the `atoms` array.
     * @param {object} atom - atom proxy
     */
    addAtom(atom) {
        if (this.indexOf(atom) < 0) this.atoms.push(atom);
        const k = atom.index;
        this.atomWeights[k] = (this.atomWeights[k] || 0) + 1;
    }

    /**
     * Decrement an atom's weight by 1, removing it from the bead entirely
     * once its weight reaches 0 (the shift-click path).
     * @param {object} atom - atom proxy
     */
    removeAtom(atom) {
        const k = atom.index;
        if (!this.atomWeights[k]) return;
        this.atomWeights[k] -= 1;
        if (this.atomWeights[k] <= 0) {
            delete this.atomWeights[k];
            const idx = this.indexOf(atom);
            if (idx >= 0) this.atoms.splice(idx, 1);
        }
    }

    /**
     * Plain-click handling: always adds/increments (see addAtom). Despite
     * the name, this never removes — shift-click decrementing goes through
     * removeAtom directly instead, so the two click modifiers map to two
     * distinct methods rather than one true toggle.
     * @param {object} atom - atom proxy
     */
    toggleAtom(atom) {
        this.addAtom(atom);
    }

    /** Display/export name. */
	set name(name) { this._name = name; }
	get name()     { return this._name; }

    /**
     * Martini type code (e.g. "SP2a"). Defaults to the placeholder "TYPe"
     * (see the constructor) — that exact string is the sentinel
     * updateMappingStats (visualization.js) checks for to decide whether
     * every bead has actually had a real type assigned yet.
     */
    set type(value) { this._type = value; }
    get type()      { return this._type; }

    /**
     * Formal charge, in units of e. The setter coerces its input with
     * `parseFloat(value) || 0`, so it accepts a string directly from an
     * HTML number input's `.value`, with an empty/non-numeric input
     * becoming 0 rather than NaN.
     */
    set charge(value) { this._charge = parseFloat(value) || 0; }
    get charge()      { return this._charge; }

    /**
     * Residue name, read from the first assigned atom (assumes a bead's
     * atoms all belong to the same residue, true for typical single-residue
     * Martini mappings). "UNK" if the bead has no atoms yet.
     */
	get resname() {
	    if (this.atoms.length < 1) return 'UNK';
	    return this.atoms[0].resname;
    }

    /** Residue number, read from the first assigned atom (see resname); 0 if empty. */
	get resid() {
	    if (this.atoms.length < 1) return 0;
	    return this.atoms[0].resno;
    }

    /**
     * Whether `atom` is currently part of this bead.
     * @param {object} atom - atom proxy
     * @returns {boolean}
     */
	isAtomIn(atom) {
		return this.indexOf(atom) >= 0;
	}

    /**
     * This bead's 3D centre: the weight-averaged position of its atoms
     * (each atom counted `atomWeights[atom.index]` times, default 1). This
     * is both where the bead's sphere is drawn and the position used by
     * sasa.js's cgSASA / fileformats.js's generateGRO / sasa.js's
     * beadsToPDB.
     * @returns {object} an NGL.Vector3
     */
    get center() {
        let mass = 0;
        let position = new NGL.Vector3(0, 0, 0);
        for (const atom of this.atoms) {
            const w = this.atomWeights[atom.index] || 1;
            mass += w;
            for (let i = 0; i < w; i++) position.add(atom.positionToVector3());
        }
        position.divideScalar(mass);
        return position;
    }

    /**
     * This bead's atoms, with each atom repeated `atomWeights[atom.index]`
     * times — a flattened view for callers that need every weighted
     * contribution rather than every distinct atom (e.g. Shaker's own
     * format expresses a weighted atom by listing its name more than once,
     * rather than via an explicit weight field — see
     * BeadCollection.expandedAtomNames).
     * @returns {Array} atom proxies, with duplicates per weight
     */
    expandedAtoms() {
        let out = [];
        for (const atom of this.atoms) {
            const w = this.atomWeights[atom.index] || 1;
            for (let i = 0; i < w; i++) out.push(atom);
        }
        return out;
    }
}


/**
 * The full bead mapping for one structure: every Bead, which one is
 * currently selected for editing, and a record of the structure's original
 * atom names (see setOriginalAtomNames).
 */
export class BeadCollection {
    constructor () {
        this._beads = [];
        this._current = null;
        this._largestIndex = -1;
        this._atomNames = new Map();
        this.newBead();
    }

    /**
     * Create a new, empty bead, auto-named "B" + a monotonic counter that
     * never decreases or gets reused — even after beads are deleted, the
     * next new bead's number keeps counting up (contrast with clearBeads,
     * which does reset the counter). Becomes the selected bead.
     * @returns {Bead} the new bead
     */
    newBead () {
        let bead = new Bead();
        this._largestIndex += 1;
        bead.name = 'B' + this._largestIndex;
        this._beads.push(bead);
        this._current = bead;
        return bead;
    }

    /**
     * Remove the bead at `index`. Does not itself adjust `currentBead` if
     * the removed bead was selected — callers are responsible for
     * re-selecting (or creating a new bead) afterward if needed.
     * @param {number} index
     */
    removeBead(index) {
        this._beads.splice(index, 1);
    }

    /** The bead currently selected for editing, or null if none is. */
    get currentBead() { return this._current; }
    /** Every bead in this collection, in creation/display order. */
    get beads()       { return this._beads; }

    /** Select the bead at `index` as the current bead for editing. */
    selectBead(index) {
        this._current = this._beads[index];
    }

    /** Clear the current bead selection (none selected). */
    deselectBead() {
        this._current = null;
    }

    /** Remove every bead and reset bead-name numbering back to B0. */
    clearBeads() {
        this._beads = [];
        this._current = null;
        this._largestIndex = -1;
    }

    /**
     * How many beads currently include `atom` — used to flag atoms shared
     * between multiple beads in the UI.
     * @param {object} atom - atom proxy
     * @returns {number}
     */
    countBeadsForAtom(atom) {
        let count = 0;
        for (const bead of this.beads) {
            if (bead.isAtomIn(atom)) count += 1;
        }
        return count;
    }

    /**
     * Record the loaded structure's original atom names, keyed by atom
     * index — populated from fileformats.js's readOriginalAtomNames, since
     * NGL's own parsed `atom.atomname` isn't guaranteed to match the
     * source file byte-for-byte. atomName falls back to NGL's name for any
     * index not present here.
     * @param {Array<string>} names - one name per atom, in atom-index order
     */
    setOriginalAtomNames(names) {
        this._atomNames.clear();
        if (!Array.isArray(names)) return;
        for (let i = 0; i < names.length; i++) {
            this._atomNames.set(i, names[i]);
        }
    }

    /**
     * The name to use for `atom` in exported output: the original
     * source-file name if recorded (see setOriginalAtomNames), otherwise
     * NGL's own parsed name.
     * @param {object} atom - atom proxy
     * @returns {string}
     */
    atomName(atom) {
        return this._atomNames.get(atom.index) ?? atom.atomname;
    }

    /**
     * atomName applied to a list of atoms.
     * @param {Array} atoms - atom proxies
     * @returns {Array<string>}
     */
    atomNames(atoms) {
        return atoms.map((atom) => this.atomName(atom));
    }

    /**
     * Names for a bead's expandedAtoms — i.e. an atom weighted ×2 has its
     * name listed twice. This is the form generatePythonAssignments
     * (fileformats.js) writes into a Shaker mapping's "atoms" list.
     * @param {Bead} bead
     * @returns {Array<string>}
     */
    expandedAtomNames(bead) {
        return bead.expandedAtoms().map((atom) => this.atomName(atom));
    }

    /**
     * atomName for every atom in a structure, in structure iteration order
     * — used for NGL label representations (see attachAALabels in
     * visualization.js) so the AA labels toggle shows the original names
     * too, not just NGL's own.
     * @param {object} structure - NGL-style structure (eachAtom)
     * @returns {Array<string>}
     */
    structureAtomNames(structure) {
        const names = [];
        if (structure && typeof structure.eachAtom === "function") {
            structure.eachAtom((atom) => names.push(this.atomName(atom)));
        }
        return names;
    }
}
