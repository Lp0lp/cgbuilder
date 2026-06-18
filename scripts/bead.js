export class Bead {
	constructor () {
		this._name = null;
        this._type = "TYPe";
		this._charge = 0;
		this.atoms = [];
        this.atomWeights = {};        // key: atom.index -> integer weight
	}

    indexOf(atom) {
        for (let i = 0; i < this.atoms.length; i++) {
            if (this.atoms[i].index === atom.index) return i;
        }
        return -1;
    }

    addAtom(atom) {
        if (this.indexOf(atom) < 0) this.atoms.push(atom);
        const k = atom.index;
        this.atomWeights[k] = (this.atomWeights[k] || 0) + 1;
    }

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

    // IMPORTANT: click -> add weight; shift-click decrements via removeAtom
    toggleAtom(atom) {
        this.addAtom(atom);
    }

	set name(name) { this._name = name; }
	get name()     { return this._name; }

    set type(value) { this._type = value; }
    get type()      { return this._type; }

    set charge(value) { this._charge = parseFloat(value) || 0; }
    get charge()      { return this._charge; }

	get resname() {
	    if (this.atoms.length < 1) return 'UNK';
	    return this.atoms[0].resname;
    }

	get resid() {
	    if (this.atoms.length < 1) return 0;
	    return this.atoms[0].resno;
    }

	isAtomIn(atom) {
		return this.indexOf(atom) >= 0;
	}

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

    expandedAtoms() {
        let out = [];
        for (const atom of this.atoms) {
            const w = this.atomWeights[atom.index] || 1;
            for (let i = 0; i < w; i++) out.push(atom);
        }
        return out;
    }
}


export class BeadCollection {
    constructor () {
        this._beads = [];
        this._current = null;
        this._largestIndex = -1;
        this._atomNames = new Map();
        this.newBead();
    }

    newBead () {
        let bead = new Bead();
        this._largestIndex += 1;
        bead.name = 'B' + this._largestIndex;
        this._beads.push(bead);
        this._current = bead;
        return bead;
    }

    removeBead(index) {
        this._beads.splice(index, 1);
    }

    get currentBead() { return this._current; }
    get beads()       { return this._beads; }

    selectBead(index) {
        this._current = this._beads[index];
    }

    clearBeads() {
        this._beads = [];
        this._current = null;
        this._largestIndex = -1;
    }

    countBeadsForAtom(atom) {
        let count = 0;
        for (const bead of this.beads) {
            if (bead.isAtomIn(atom)) count += 1;
        }
        return count;
    }

    setOriginalAtomNames(names) {
        this._atomNames.clear();
        if (!Array.isArray(names)) return;
        for (let i = 0; i < names.length; i++) {
            this._atomNames.set(i, names[i]);
        }
    }

    atomName(atom) {
        return this._atomNames.get(atom.index) ?? atom.atomname;
    }

    atomNames(atoms) {
        return atoms.map((atom) => this.atomName(atom));
    }

    expandedAtomNames(bead) {
        const names = [];
        for (const atom of bead.atoms) {
            const weight = bead.atomWeights[atom.index] || 1;
            for (let i = 0; i < weight; i++) names.push(this.atomName(atom));
        }
        return names;
    }

    structureAtomNames(structure) {
        const names = [];
        if (structure && typeof structure.eachAtom === "function") {
            structure.eachAtom((atom) => names.push(this.atomName(atom)));
        }
        return names;
    }
}
