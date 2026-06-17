class Bead {
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
        // first time: add to unique list
        if (this.indexOf(atom) < 0) this.atoms.push(atom);
        // always increment weight
        const k = atom.index;
        this.atomWeights[k] = (this.atomWeights[k] || 0) + 1;
    }

    removeAtom(atom) {
        const k = atom.index;
        if (!this.atomWeights[k]) return;

        this.atomWeights[k] -= 1;

        // if weight reaches 0, remove atom from unique list
        if (this.atomWeights[k] <= 0) {
        delete this.atomWeights[k];
        const idx = this.indexOf(atom);
        if (idx >= 0) this.atoms.splice(idx, 1);
        }
    }

	set name(name) {
		this._name = name;
	}

	get name() {
		return this._name;
	}

	get resname() {
	    if (this.atoms.length < 1) {
	        return 'UNK';
        }
	    return this.atoms[0].resname;
    }

	get resid() {
	    if (this.atoms.length < 1) {
	        return 0;
        }
	    return this.atoms[0].resno;
    }

    set type(value) {
        this._type = value;
    }

    get type() {
        return this._type;
    }

    set charge(value) {
        this._charge = parseFloat(value) || 0;
    }

    get charge() {
        return this._charge;
    }

	isAtomIn(atom) {
		return this.indexOf(atom) >= 0;
	}

    // IMPORTANT: change toggle behavior
    // click -> add weight; (you can add a separate "decrement" action)
    toggleAtom(atom) {
        this.addAtom(atom);
    }

    // weighted center
    get center() {
        let mass = 0;
        let position = new NGL.Vector3(0, 0, 0);
        for (const atom of this.atoms) {
        const w = this.atomWeights[atom.index] || 1;
        mass += w;
        // add atom position w times (cheap way without needing vector scaling)
        for (let i = 0; i < w; i++) position.add(atom.positionToVector3());
        }
        position.divideScalar(mass);
        return position;
    }

    // helper to export duplicates (for your python output)
    expandedAtoms() {
        let out = [];
        for (const atom of this.atoms) {
        const w = this.atomWeights[atom.index] || 1;
        for (let i = 0; i < w; i++) out.push(atom);
        }
        return out;
    }
}


class BeadCollection {
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

    get currentBead() {
        return this._current;
    }

    get beads() {
        return this._beads;
    }

    selectBead(index) {
        this._current = this._beads[index];
    }

    countBeadsForAtom(atom) {
        let count = 0;
        for (const bead of this.beads) {
            if (bead.isAtomIn(atom)) {
                count += 1;
            }
        }
        return count;
    }

    setOriginalAtomNames(names) {
        this._atomNames.clear();
        if (!Array.isArray(names)) {
            return;
        }
        // Parsed input atom names must remain aligned with NGL atom.index.
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
            for (let i = 0; i < weight; i++) {
                names.push(this.atomName(atom));
            }
        }
        return names;
    }

    structureAtomNames(structure) {
        const names = [];
        if (structure && typeof structure.eachAtom === "function") {
            structure.eachAtom((atom) => {
                names.push(this.atomName(atom));
            });
        }
        return names;
    }
}


class Visualization {
    constructor(collection, stage) {
        this.collection = collection;
        this.representation = null;
        this.stage = stage;
        this.shapeComp = null;
        this.showCG = false;

        // Solvent-accessible surface area (SASA) state.
        this.component = null;        // the loaded AA structure component
        this.aaSurface = null;        // NGL surface representation on the AA structure
        this.showAASurface = false;
        this.cgSurfaceComp = null;    // synthetic bead component carrying the CG surface
        this.showCGSurface = false;
        this._cgSurfaceToken = 0;     // guards against stale async surface loads
        this._aaSASAValue = null;     // cached; recomputed only on molecule load

        let toggleCG = document.getElementById('toggle-cg');
        toggleCG.onclick = (event) => this.onToggleCG(event);
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

        document.getElementById('dl-ndx').onclick = (event) => {
            download('cgbuilder.ndx', generateNDX(this.collection));
        };
        document.getElementById('dl-map').onclick = (event) => {
            download('cgbuilder.map', generateMap(this.collection));
        };
        document.getElementById('dl-gro').onclick = (event) => {
            download('cgbuilder.gro', generateGRO(this.collection));
        };
        document.getElementById('dl-py').onclick = (event) => {
            download('cgbuilder_assignments.py', generatePythonAssignments(this.collection));
        };
        document.getElementById('copy-py').onclick = async (event) => {

            const button = event.target;
            const originalText = button.textContent;

            const text = document.getElementById('py-output').textContent || "";

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
                setTimeout(() => {
                    button.textContent = originalText;
                }, 1200);
            }
        };
    }

	get currentBead() {
	    return this.collection.currentBead;
	}

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

        // Solvent-accessible surface of the all-atom structure. NGL derives the
        // per-atom vdW radii from its built-in element table.
        this.aaSurface = component.addRepresentation("surface", {
            surfaceType: "sas",
            probeRadius: PROBE_RADIUS,
            color: "#f4b642",
            opacity: 0.3,
            visible: this.showAASurface,
            useWorker: false,
        });

        // Pre-compute AA SASA once; the structure doesn't change after load.
        this._aaSASAValue = aaSASA(component.structure, PROBE_RADIUS);
        this.updateSASA();
    }

    attachAALabels(component) {
        this.aa_labels = component.addRepresentation(
            "label",
            {
                labelType: "text",
                labelText: this.collection.structureAtomNames(component.structure),
                labelGrouping: "atom",
            },
        );

        let buttons = document.getElementsByClassName("toggle-aa-labels");
        for (const button of buttons) {
            button.disabled = false;
            button.onclick = (event) => this.onToggleAALabels(event);
        }
    }

    onToggleCG(event) {
        this.showCG = (! this.showCG);
        this.drawCG();
    }

    onToggleAASurface(event) {
        this.showAASurface = event.target.checked;
        if (this.aaSurface) {
            this.aaSurface.setVisibility(this.showAASurface);
        }
    }

    onToggleCGSurface(event) {
        this.showCGSurface = event.target.checked;
        this.drawCGSurface();
    }

    async onPredictTypes() {
        const btn = document.getElementById('predict-types');
        const originalText = btn.textContent;
        btn.disabled = true;

        try {
            btn.textContent = 'Loading RDKit…';
            const RDKit = await loadRDKit();

            // Pre-canonicalize the FRAG_DELTA_F table keys once so lookups work
            // regardless of whether the fragment SMILES is aromatic or Kekulé.
            const canonTable = buildCanonTable(RDKit);

            btn.textContent = 'Predicting…';
            const hasExplicitH = structureHasHydrogens(
                this.component && this.component.structure);
            // Whole-molecule chemistry perception from geometry: aromatic rings
            // and multiple-bond orders, so carbonyls/amides/aromatics are read
            // correctly even from PDB/GRO files that carry no bond-order data.
            const { aromaticAtoms, bondOrders } = perceiveChemistry(
                this.component && this.component.structure);

            for (const bead of this.collection.beads) {
                if (bead.atoms.length === 0) continue;

                const smiles = fragmentToSmiles(bead.atoms, { aromaticAtoms, bondOrders });
                if (!smiles) continue;

                const mol = RDKit.get_mol(smiles);
                if (!mol) {
                    console.warn(`Bead "${bead.name}": invalid SMILES "${smiles}" — skipping`);
                    continue;
                }

                // Canonicalize so it matches the pre-canonicalized table keys.
                const canonSmiles = mol.get_smiles();
                const desc = JSON.parse(mol.get_descriptors());
                mol.delete();

                const charge     = bead.atoms.reduce((s, a) => s + (a.formalCharge ?? 0), 0)
                                   + (bead.charge || 0);
                const hasHalogen = bead.atoms.some(a =>
                    ['F','CL','BR','I'].includes((a.element || '').toUpperCase()));
                const inRing     = bead.atoms.some(a => aromaticAtoms.has(a.index));
                const heavyCount = bead.atoms.filter(a =>
                    (a.element || 'C').toUpperCase() !== 'H').length;
                const hDonors    = beadDonorCount(bead, hasExplicitH);

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
                        if (!aromaticAtoms.has(a.index)) continue;
                        const k = fragmentToSmiles(bead.atoms, {
                            aromaticNotation: true, aromaticAtoms, bondOrders,
                            startIndex: a.index });
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
                    // RDKit.js MinimalLib key is CrippenClogP (not MolLogP).
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
                    { deltaF, charge, hDonors, hAcceptors, hasHalogen, inRing, heavyCount });
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

    // Build (or refresh) the CG solvent-accessible surface. Each bead is written
    // to a tiny synthetic PDB with its vdW radius stored in the B-factor column,
    // then NGL's surface generator reads that radius via radiusType: "bfactor".
    drawCGSurface() {
        if (this.cgSurfaceComp != null) {
            this.stage.removeComponent(this.cgSurfaceComp);
            this.cgSurfaceComp = null;
        }
        if (! this.showCGSurface) {
            return;
        }
        let pdb = beadsToPDB(this.collection);
        if (! pdb) {
            return;
        }
        let token = ++this._cgSurfaceToken;
        this.stage
            .loadFile(new Blob([pdb], {type: "text/plain"}), {ext: "pdb"})
            .then((comp) => {
                // Drop the result if a newer request superseded this one or the
                // surface was toggled off while loading.
                if (token !== this._cgSurfaceToken || ! this.showCGSurface) {
                    this.stage.removeComponent(comp);
                    return;
                }
                comp.addRepresentation("surface", {
                    surfaceType: "sas",
                    radiusType: "bfactor",
                    radiusScale: 1.0,
                    probeRadius: PROBE_RADIUS,
                    color: "#7fc8a9",
                    opacity: 0.3,
                    useWorker: false,
                });
                this.cgSurfaceComp = comp;
            })
            .catch((err) => {
                console.error("Error building CG surface:", err);
            });
    }

    onToggleAALabels(event) {
        let visible = ! this.aa_labels.visible;
        this.aa_labels.setVisibility(visible);
        let text;
        if (visible) {
            text = 'Hide labels';
        } else {
            text = 'Show labels';
        }
        let buttons = document.getElementsByClassName("toggle-aa-labels");
        for (const button of buttons) {
            button.textContent = text;
        }
    }

    onClick(pickingProxy) {
    if (pickingProxy && pickingProxy.atom) {
        if (pickingProxy.mouse && pickingProxy.mouse.shiftKey) {
        this.currentBead.removeAtom(pickingProxy.atom);  // decrement
        } else {
        this.currentBead.addAtom(pickingProxy.atom);     // increment
        }
        this.updateSelection();
    }
    }

	onNewBead(event) {
	    this.collection.newBead();
	    this.updateSelection();
	}

    onBeadSelected(event) {

        const tag = event.target.tagName;

        if (tag === "INPUT" || tag === "BUTTON" || tag === "FORM" || tag === "LABEL") {
            return;
        }

        let realTarget = findParentWithClass(event.target, "bead-view");
        let nodes = document.getElementById("bead-list").childNodes;
        let index = 0;

        for (const child of nodes) {
            if (child === realTarget) {
                this.collection.selectBead(index);
            }
            index += 1;
        }

        this.updateSelection();
    }

	onBeadRemove(event) {
        let realTarget = findParentWithClass(event.target, "bead-view");
        let nodes = document.getElementById("bead-list").childNodes;
        let index = 0;
        let selected = -1;
        for (const child of nodes) {
            if (child === realTarget) {
                selected = index;
                break;
            }
            index += 1;
        }
        if (selected >= 0) {
            this.collection.removeBead(selected);
            if (this.collection.beads.length === 0) {
                this.collection.newBead();
            }
            if (realTarget.classList.contains('selected-bead')) {
                this.collection.selectBead(0);
            }
        }

        this.updateSelection();
    }

    onNameChange(event) {
        let realTarget = findParentWithClass(event.target, "bead-view");
        let nodes = document.getElementById("bead-list").childNodes;
        let index = 0;
        for (const child of nodes) {
            if (child === realTarget) {
                this.collection.beads[index].name = event.target.value;
            }
            index += 1;
        }
        this.updateName();
    }

	selectionString(bead) {
        if (bead.atoms.length > 0) {
            let sel = "@";
            for (let i=0; i < bead.atoms.length; i++) {
                if (sel !== '@') {
                    sel = sel + ',';
                }
                sel = sel + bead.atoms[i].index;
            }
            return sel;
        }
        return "not all";
    }

    updateName() {
        this.updateNDX();
        this.updateMap();
        this.updateGRO();
        this.updatePY();
        this.updateSASA();
        this.drawCG();
    }

    updateSelection() {
        let selString = this.selectionString(this.currentBead);
        this.representation.setSelection(selString);
        this.clearBeadList();
        this.createBeadList();
        this.updateName();
    }

    createBeadListItem(bead) {

        let textNode;
        let list = document.getElementById("bead-list");
        let item = document.createElement("li");

        item.classList.add("bead-view");

       /* ===============================
        HEADER ROW (Fields + Delete)
        =============================== */

        let headerRow = document.createElement("div");
        headerRow.classList.add("bead-header");

        let fieldsNode = document.createElement("div");
        fieldsNode.classList.add("bead-fields");

        // helper to build labeled field
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
        nameNode.oninput = (event) => this.onNameChange(event);
        nameNode.addEventListener("mousedown", e => e.stopPropagation());
        addLabeledField("Name", nameNode);

        // TYPE
        let typeNode = document.createElement("input");
        typeNode.type = "text";
        typeNode.value = bead.type;
        typeNode.classList.add("bead-type");
        typeNode.oninput = (event) => {
            bead.type = event.target.value;
            this.updateName();
        };
        typeNode.addEventListener("mousedown", e => e.stopPropagation());

        // Suggestion chip — shown when a prediction has been run
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
        chargeNode.oninput = (event) => {
            bead.charge = event.target.value;
            this.updateName();
        };
        chargeNode.addEventListener("mousedown", e => e.stopPropagation());
        addLabeledField("Charge", chargeNode);

        // DELETE BUTTON
        let removeNode = document.createElement("button");
        removeNode.textContent = "Delete";
        removeNode.classList.add("delete-bead");
        removeNode.onclick = (event) => {
            event.stopPropagation();
            this.onBeadRemove(event);
        };

        // Assemble header row
        headerRow.appendChild(fieldsNode);
        headerRow.appendChild(removeNode);

        item.appendChild(headerRow);

        /* ===============================
        ATOM LIST
        =============================== */

        let nameList = document.createElement("ul");

        if (bead.atoms.length > 0) {

            for (let i = 0; i < bead.atoms.length; i++) {

                const atom = bead.atoms[i];
                const name = this.collection.atomName(atom);

                const w = (bead.atomWeights && bead.atomWeights[atom.index])
                    ? bead.atomWeights[atom.index]
                    : 1;

                let subitem = document.createElement("li");

                const label = (w > 1) ? `${name} ×${w}` : name;
                textNode = document.createTextNode(label);

                subitem.appendChild(textNode);

                if (this.collection.countBeadsForAtom(atom) > 1) {
                    let shareitem = document.createElement("abbr");
                    shareitem.title = "This atom is shared between multiple beads.";
                    shareitem.textContent = " 🔗";
                    subitem.appendChild(shareitem);
                }

                nameList.appendChild(subitem);
            }
        }

        item.appendChild(nameList);

        /* ===============================
        SELECTION HANDLING
        =============================== */

        item.onclick = (event) => this.onBeadSelected(event);

        list.appendChild(item);

        if (bead === this.currentBead) {
            item.classList.add("selected-bead");
            item.scrollIntoView({ block: 'nearest' });
        }
    }

    createBeadList() {
        for (let bead of this.collection.beads) {
            this.createBeadListItem(bead);
        }
    }

    clearBeadList() {
        let list = document.getElementById('bead-list');
        while (list.lastChild) {
            list.removeChild(list.lastChild);
        }
    }

    updateNDX() {
        let displayNode = document.getElementById('ndx-output');
        displayNode.textContent = generateNDX(this.collection);
    }

    updateMap() {
        let displayNode = document.getElementById('map-output');
        displayNode.textContent = generateMap(this.collection);
    }

    updateGRO() {
        let displayNode = document.getElementById('gro-output');
        displayNode.textContent = generateGRO(this.collection);
    }

    updatePY() {
        let displayNode = document.getElementById('py-output');
        displayNode.textContent = generatePythonAssignments(this.collection);
    }

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

    drawCG() {
        let normalColor = [0.58, 0.79, 0.66];
        let selectedColor = [0.25, 0.84, 0.96];
        let color = normalColor;
        let opacity = 0.2;
        if (this.showCG) {
            opacity = 1;
        }
        if (this.shapeComp != null) {
            this.stage.removeComponent(this.shapeComp);
        }
        let shape = new NGL.Shape("shape");
        for (let bead of this.collection.beads) {
            color = normalColor;
            if (bead === this.currentBead) {
                color = selectedColor;
            }
            if (bead.atoms.length > 0) {
                shape.addSphere(bead.center, color, 1.12, bead.name);
            }
        }
        this.shapeComp = this.stage.addComponentFromObject(shape);
        this.shapeComp.addRepresentation("buffer", {opacity: opacity});

        // Keep the CG surface in sync with the current mapping when it is shown.
        this.drawCGSurface();
    }
}


function findParentWithClass(element, className) {
    let node = element;
    while (node) {
        if (node.classList.contains(className)) {
            return node;
        }
        node = node.parentElement;
    }
    return null;
}


function generateNDX(collection) {
    let ndx = "";
    for (const bead of collection.beads) {
        ndx += "[ " + bead.name + " ]\n";
        for (const atom of bead.atoms) {
            ndx += (atom.index + 1) + " ";
        }
        ndx += "\n\n";
    }
    return ndx;
}


function generateMap(collection) {
    let output = "[ to ]\nmartini\n\n[ martini ]\n";
    let atomToBeads = {};
    let atoms = [];
    let atomname;
    let index;
    for (const bead of collection.beads) {
        output += bead.name + " ";
        for (const atom of bead.atoms) {
            atomname = collection.atomName(atom);
            if (!atomToBeads[atomname]) {
                atomToBeads[atomname] = [];
                atoms.push(atom);
            }
            atomToBeads[atomname].push(bead.name);
        }
    }
    output += "\n\n";

    output += "[ atoms ]\n";
    index = 0;
    atoms.sort(function(a, b) {return a.index - b.index});
    for (const atom of atoms) {
        index += 1;
        atomname = collection.atomName(atom);
        output += index + "\t" + atomname;
        for (const bead of atomToBeads[atomname]) {
            output += "\t" + bead;
        }
        output += "\n";
    }

    return output;
}

function generatePythonAssignments(collection) {

    const beads = collection.beads || [];
    if (beads.length === 0) return "";

    // Determine residue name
    let resname = "UNK";
    for (const bead of beads) {
        if (bead.atoms && bead.atoms.length > 0) {
            resname = bead.resname;
            break;
        }
    }

    let lines = [];

    lines.push("mapping = {");
    lines.push("    ## resname");
    lines.push(`    "${resname}": {`);
    lines.push("    ## bead name; type(opt.);   charge(opt.);         Mapping.");

    for (const bead of beads) {

        const beadName = bead.name;
        const beadType = bead.type || "type";
        const beadCharge = bead.charge ?? 0;

        const atomNames = collection.expandedAtomNames(bead)
            .map((name) => `'${name}'`)
            .join(",");

        lines.push(
            `        "${beadName}": {"type": "${beadType}", "charge": ${beadCharge}, "atoms": [${atomNames}]},`
        );
    }

    lines.push("    },");
    lines.push("}");

    return lines.join("\n") + "\n";
}

function generateGRO(collection) {
    let resid = "    0";
    let resname = "";
    let atomname = "    0";
    let atomid = 0;
    let x;
    let y;
    let z;
    let center;
    let output = "Generated with cgbuilder\n" + collection.beads.length + "\n";
    let counter = 0;
    for (const bead of collection.beads) {
        counter += 1;
        resid = String(bead.resid).padStart(5);
        atomid = String(counter).padStart(5);
        resname = bead.resname.padEnd(5).substring(0, 5);
        atomname = bead.name.padStart(5).substring(0, 5);
        center = bead.center;
        x = (center.x / 10).toFixed(3).padStart(8);
        y = (center.y / 10).toFixed(3).padStart(8);
        z = (center.z / 10).toFixed(3).padStart(8);
        output += resid + resname + atomname + atomid + x + y + z + '\n';
    }
    output += "10 10 10\n";
    return output;
}

/* Taken from <https://ourcodeworld.com/articles/read/189/how-to-create-a-file-and-generate-a-download-with-javascript-in-the-browser-without-a-server> */
function download(filename, text) {
  let element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);

  element.click();

  document.body.removeChild(element);
}

function copyTextToClipboard(text) {
    // Modern clipboard API
    if (navigator.clipboard && window.isSecureContext) {
        return navigator.clipboard.writeText(text);
    }

    // Fallback for non-HTTPS / older browsers
    let textarea = document.createElement("textarea");
    textarea.value = text;
    textarea.style.position = "fixed";
    textarea.style.left = "-9999px";
    textarea.style.top = "-9999px";
    document.body.appendChild(textarea);
    textarea.focus();
    textarea.select();
    try {
        document.execCommand("copy");
    } finally {
        document.body.removeChild(textarea);
    }
    return Promise.resolve();
}

function parseOriginalAtomNames(content, filename) {
    const lower = (filename || "").toLowerCase();
    if (lower.endsWith(".pdb") || lower.endsWith(".ent")) {
        return parsePDBAtomNames(content);
    }
    if (lower.endsWith(".gro")) {
        return parseGROAtomNames(content);
    }
    return [];
}

function parsePDBAtomNames(content) {
    const names = [];
    const lines = content.split(/\r?\n/);
    for (const line of lines) {
        if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
            names.push(line.substring(12, 16).trim());
        }
    }
    return names;
}

function parseGROAtomNames(content) {
    const names = [];
    const lines = content.split(/\r?\n/);
    if (lines.length < 3) {
        return names;
    }
    const count = parseInt(lines[1].trim(), 10);
    if (!Number.isFinite(count) || count <= 0) {
        return names;
    }
    const atomLines = lines.slice(2, 2 + count);
    for (const line of atomLines) {
        if (line.length >= 15) {
            names.push(line.substring(10, 15).trim());
        }
    }
    return names;
}

function readOriginalAtomNames(file) {
    const lower = (file && file.name ? file.name : "").toLowerCase();
    // Skip compressed input: NGL can read it, but browser text parsing here cannot.
    if (lower.endsWith(".gz")) {
        return Promise.resolve([]);
    }
    return file.text()
        .then((content) => parseOriginalAtomNames(content, file.name))
        .catch(() => []);
}

// Helper to ensure that double and triple bonds are properly rendered
function bondAwareRepresentationParams(overrides = {}) {
    return Object.assign(
        {
            multipleBond: true,
            bondSpacing: 1,
            bondScale: 0.4,
        },
        overrides
    );
}

/* ===========================================================================
   Solvent-accessible surface area (SASA)
   ===========================================================================
   Probe radius in Angstrom. 0.191 nm = 1.91 Å is the Martini tiny-bead radius,
   used here so the AA and CG surfaces are computed with the same probe and are
   directly comparable. */
const PROBE_RADIUS = 1.91;

/* Martini bead vdW radii in Angstrom, keyed by size class. Values come from the
   standard Martini bead radii (nm -> A): regular 0.264, small 0.230, tiny 0.191.
   "U" beads are virtual/ghost beads and contribute no surface (radius 0).
   Edit this table to switch force fields. */
const BEAD_RADII = {
    R: 2.64,
    S: 2.30,
    T: 1.91,
    U: 0.0,
};
const DEFAULT_BEAD_SIZE = "R"; // fall back to a regular bead for unknown types

// Resolve the Martini size class from a bead type string. In Martini naming the
// leading letter encodes the size: "S" (small) and "T" (tiny) prefix the
// chemical class (e.g. "SP2", "TC3"); "U" marks a virtual bead; everything else
// is a regular bead.
function beadSizeClass(type) {
    if (!type) {
        return DEFAULT_BEAD_SIZE;
    }
    const first = type.trim().charAt(0).toUpperCase();
    if (first === "S" || first === "T" || first === "U") {
        return first;
    }
    return DEFAULT_BEAD_SIZE;
}

function beadRadius(bead) {
    return BEAD_RADII[beadSizeClass(bead.type)];
}

// Pad an atom name into PDB columns 13-16. Names shorter than 4 characters get a
// leading space, matching the convention NGL expects when parsing.
function formatPDBAtomName(name) {
    name = (name || "").substring(0, 4);
    if (name.length >= 4) {
        return name;
    }
    return (" " + name).padEnd(4);
}

// Serialise the CG beads to a minimal PDB string, one ATOM per bead, with the
// bead vdW radius stored in the temperature-factor column (cols 61-66). Beads
// with no atoms or a zero radius (virtual "U" beads) are skipped.
function beadsToPDB(collection) {
    let lines = [];
    let serial = 0;
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) {
            continue;
        }
        const radius = beadRadius(bead);
        if (radius <= 0) {
            continue;
        }
        serial += 1;
        const center = bead.center;
        const serStr = String(serial % 100000).padStart(5);
        const name = formatPDBAtomName(bead.name);
        const resname = (bead.resname || "BEA").substring(0, 3).padEnd(3);
        const resid = String(((bead.resid % 10000) + 10000) % 10000).padStart(4);
        const x = center.x.toFixed(3).padStart(8);
        const y = center.y.toFixed(3).padStart(8);
        const z = center.z.toFixed(3).padStart(8);
        const occ = "  1.00";
        const bfac = radius.toFixed(2).padStart(6);
        lines.push(
            "ATOM  " +   // 1-6   record name
            serStr +     // 7-11  serial
            " " +        // 12    (blank)
            name +       // 13-16 atom name
            " " +        // 17    altLoc
            resname +    // 18-20 resName
            " " +        // 21    (blank)
            "A" +        // 22    chainID
            resid +      // 23-26 resSeq
            " " +        // 27    iCode
            "   " +      // 28-30 (blank)
            x + y + z +  // 31-54 coordinates
            occ +        // 55-60 occupancy
            bfac +       // 61-66 tempFactor (used as radius)
            "          " + // 67-76 (blank)
            " C"         // 77-78 element
        );
    }
    if (lines.length === 0) {
        return "";
    }
    lines.push("END");
    return lines.join("\n") + "\n";
}

/* ===========================================================================
   Bead type prediction — AutoMartini M3 port
   ===========================================================================
   Fragment → delta_f (kJ/mol) lookup table — the full AutoMartini M3
   logP_smi_extended.dat (data/logP_smi_extended.dat), with duplicate keys
   resolved by FIRST occurrence (the original calibrated value; later duplicates
   are AutoMartini's dict-overwrite bug). Keys are kept verbatim in AutoMartini's
   notation: lowercase = aromatic atoms, uppercase = aliphatic; many are
   open-chain aromatic SMILES ("cc", "cn", "ncs") that RDKit cannot parse — these
   are matched directly against the aromatic-notation fragment key (see
   buildCanonTable / onPredictTypes). Do NOT "canonicalise" or hand-edit values. */
const FRAG_DELTA_F = {
    "CC":12.0, "CCC":14.2, "CCCC":18.9, "CC(C)(C)C":18.9,
    "C=C":6.4, "C=CC":8.4, "CC=C":8.4, "CCc":8.4,
    "cCC":8.4, "cc(c)C":13.4, "cc(c)-c":13.4, "CCCc":13.4,
    "CC(c)C":13.4, "ccCC":13.4, "C=CCC":13.4, "CC=CC":13.4,
    "CCC=C":13.4, "CC1(C)CC1":13.4, "cC1CC1":13.8, "cC(C)C":13.8,
    "C=C(C)C":13.8, "C=C(C)c":13.8, "CC(C)c":13.8, "CC(C)(C)c":13.8,
    "CC(C)C":14.8, "CCC(C)C":14.8, "cC":6.4, "Cc":6.4,
    "ccC":8.4, "cc":4.5, "Ccc":6.3,
    "ccc":6.3, "cc-c":6.3, "C=Cc":6.3, "C1CC1":9.2,
    "cCcc":11.2, "CC1CC1":11.2, "Cc(c)c":11.2, "cccc":11.2,
    "cccC":11.2, "C#C":3.6, "C#CC":5.3, "C#Cc":5.3,
    "CC#C":5.3, "C#CCC":10.1, "CC#CC":10.1, "CCC#C":10.1,
    "C#CC(C)C":10.1, "CS":3.6, "CCS":5.3, "CCCS":10.1,
    "cs":3.6, "cS":3.6, "ccs":5.3, "ccS":5.3,
    "cSc":5.3, "Scc":5.3, "Ccs":5.3, "c(CS)c":10.1,
    "C=C(C)S":10.1, "SC":3.6, "CSC":5.3, "SCC":5.3,
    "csc":5.3, "C=CS":5.3, "CSc":8.4, "CSCC":10.1,
    "cc(-c)s":10.1, "CCSC":10.1, "SCCC":10.1, "cc(C)s":10.1,
    "CC(C)S":13.4, "CC(C)(C)S":13.4, "CCF":2.7, "CCCF":4.3,
    "CCCCF":8.7, "FCC":2.7, "C(F)CC":4.3, "CC(F)C":4.3,
    "CC(F)CC":8.7, "cCl":5.4, "ccCl":8.0, "scCl":8.0,
    "CCl":5.4, "CCCl":8.0, "CCCCl":13.9, "ClCC":8.0,
    "C=CCl":8.0, "SCCl":8.0, "ocCl":4.3, "C=CF":9.4,
    "ncF":9.4, "C=C(c)Cl":13.9, "CC(Cl)C":13.9, "CC(C)Cl":13.9,
    "cc(c)Cl":13.9, "cc(s)Cl":13.9, "CC(F)F":13.9, "CC(C)(C)Cl":14.3,
    "CBr":5.2, "CCBr":7.2, "CCCBr":12.7, "cBr":5.2,
    "BrCC":7.2, "C=CBr":7.2, "ccBr":7.2, "ClCCl":7.2,
    "CC(Br)C":12.7, "C=C(C)Br":12.7, "C=C(C)Cl":12.7, "NC(F)F":12.7,
    "FC(F)S":12.7, "CC(C)(F)F":13.9, "cc(c)F":14.3, "CC(C)F":14.3,
    "OC(F)F":8.7, "CC(F)(F)F":8.7, "ClC1(Cl)CC1":8.7, "CC(Cl)Cl":8.7,
    "cC(F)(F)F":8.7, "OC(F)(F)F":8.7, "CC(O)F":8.7, "OCCl":4.3,
    "O=PCl":4.3, "FC(F)F":8.7, "I":7.6, "CF":7.6,
    "cF":7.6, "CI":7.6, "CCI":14.3, "IC":7.6,
    "cI":7.6, "ccI":9.4, "FCF":7.2, "ccF":8.0,
    "C(I)C":14.3, "Cn":2.3, "cN":-2.8, "CN":-4.1,
    "CCN":-2.5, "CC#N":0.2, "cC#N":1.1, "C=CN":-2.5,
    "CNC":3.6, "CCn":3.6, "CCCN":1.0, "CN(C)C":8.1,
    "CN(C)c":8.1, "CC(C)n":8.1, "CC(c)N":8.1, "CCC(C)(C)N":10.7,
    "NS":-4.1, "cn":-4.1, "C=N":-4.1,
    "O=s":-4.1, "NC=O":-6.9, "NCO":-5.4, "NCC":-0.9,
    "Ccn":6.0, "C#N":-1.2, "ccN":-2.5, "C=Nc":-2.5,
    "cNc":-2.5, "ccn":-1.0, "C=Cn":-1.0, "cnc":-1.0,
    "ncs":-1.0, "cns":-1.0, "NCS":-1.0, "cSC":5.3,
    "cCS":5.3, "C=NC":3.6, "cNC":3.6, "cno":1.1,
    "C=NO":1.1, "nco":1.1, "Nco":1.1, "cC=O":1.1,
    "con":1.1, "cOc":1.1, "nCO":-5.2, "Cnn":-5.2,
    "CNN":-5.2, "cNN":-5.2, "Ncs":-2.5, "O=CS":-2.5,
    "Cs=O":-2.5, "O=s=O":-3.7, "ncn":-3.7, "cnn":-3.7,
    "c-nn":-3.7, "nns":-3.7, "ncN":-6.9, "cnN":-6.9,
    "N=CN":-6.9, "C=Nn":-6.9, "CNc":-6.9, "c-cn":-2.5,
    "cnC":-2.5, "cn-c":-2.5, "CC=N":-2.5, "N=CS":-2.5,
    "ncS":-2.5, "csn":-2.5, "COn":-2.5, "C=NN":-1.0,
    "NCN":-5.4, "cCN":-4.2, "CNO":0.2, "nn":-6.1,
    "nN":-6.1, "Nn":-6.1, "NN":-7.8, "ns":-4.1,
    "no":-7.2, "N=O":-7.2, "O=P":-5.0, "NO":-5.0,
    "N=S=O":-1.0, "ccnn":0.2, "cncn":0.2, "nc(n)S":0.2,
    "Cc(n)n":0.2, "cc(n)n":0.2, "cn(-c)n":0.2, "nc(n)N":0.2,
    "CN(c)N":0.2, "C=C(N)N":0.2, "N[SH](=O)=O":-1.9, "CS(N)(=O)=O":-1.9,
    "N=C(N)N":-1.9, "c[N+](=O)[O-]":3.8, "CC(C)nn":3.8, "C[SH](=O)=O":0.2,
    "c[SH](=O)=O":0.2, "CS(C)(=O)=O":0.2, "N=[SH](N)=O":0.2, "NC(=O)S":0.2,
    "Cc(n)s":0.2, "c-c(n)s":0.2, "cc(N)s":0.2, "CS(C)=O":0.2,
    "cc(c)n":1.0, "CC(=N)S":1.0, "cc(n)S":1.0, "cc(-c)n":1.0,
    "C=C(C)n":1.0, "cn(c)C":1.0, "C#CCN":1.0, "CC(C)(C)n":8.1,
    "nc(N)s":-1.9, "cc(n)N":-1.9, "CC(=N)N":-1.9, "cn(C)n":-1.9,
    "N=C(N)S":-1.9, "c-c(n)n":-1.9, "ccC#N":4.3, "CC(C)C#N":4.3,
    "NC(N)=O":0.2, "C=C(N)O":0.2, "ncO":-5.2, "N=NN":-5.2,
    "CON":-1.0, "NC(=O)":-3.7, "N=CO":-3.7, "N[SH]=O":-3.7,
    "O=[N+][O-]":-1.8, "CC(C)=N":2.2, "cc(C)n":2.2, "CC(C)N":1.0,
    "CNC(C)C":1.0, "cN(C)C":1.0, "cC(C)N":1.0, "cccn":1.0,
    "ccnc":1.0, "nc(N)o":1.0, "CC(N)C":2.2, "C=C(C)N":2.2,
    "CC(N)=O":1.0, "CC(=N)O":1.0, "cC(C)=N":1.0, "cc(n)s":0.2,
    "cC(N)=O":0.2, "nc(n)O":0.2, "N=C(N)O":0.2, "C=C(N)n":0.2,
    "Cc(n)o":4.3, "cc(n)O":-1.9, "c-c(n)o":-1.9, "ccon":-1.1,
    "nC1CC1":8.1, "NC1CC1":8.1, "CCNC":8.1, "CCC(C)N":10.7,
    "cc(c)S":10.1, "cc(-c)o":7.8, "NCCO":7.8, "CC(C)(C)O":8.1,
    "CC(C)(C)N":8.1, "ccN(C)C":8.1, "CC(C)(c)O":8.1, "cc(c)N":3.8,
    "CCC#N":3.8, "CN(c)C":2.2, "ccNC":2.2, "CC(C)(N)O":4.3,
    "OC":-1.4, "CO":-1.4, "co":2.3, "COc":3.8,
    "COC":0.2, "CCO":0.2, "Cco":0.2, "cCO":-4.2,
    "cco":3.8, "CCOC":3.8, "COCC":3.8, "coc":3.8,
    "CP=O":-2.5, "C=O":-2.8, "CC=O":-1.0, "cC=N":-1.0,
    "CCC=O":2.2, "cC(O)O":2.2, "CC(C)O":9.9, "CC(c)O":-2.0,
    "cC(C)O":-2.0, "OCC":-1.0, "O=S=O":-3.7, "CC(=O)C":2.2,
    "C(=O)O":-1.2, "O=CO":-6.9, "CC(O)O":1.0, "CC(=O)O":-3.8,
    "cC(=O)O":-3.8, "ccOC":7.8, "COC(C)C":7.8, "cc(C)o":7.8,
    "C=C(C)O":7.8, "cc(c)o":7.8, "ccOc":7.8, "ccCO":4.3,
    "cOC":1.1, "CCCO":-2.0, "CC(O)C":-2.0, "cO":-6.1,
    "Ocs":-4.2, "ccO":-4.2, "C=CO":-4.2, "cc(c)O":-4.2,
    "cccO":-4.2,
};

/* Free energies of transfer (kJ/mol) for all Martini 3 bead types.
   Source: Souza et al. 2021 (https://doi.org/10.1038/s41592-021-01098-3). */
const DELTA_F = {
    C1:18.9, C2:14.8, C3:13.8, C4:13.4, C5:11.2, C6:10.1,
    N1:8.1,  N2:5.6,  N3:1.8,  N4:2.2,  N5:0.0,  N6:-1.1,
    P1:-2.0, P2:-3.8, P3:-5.1, P4:-7.4, P5:-9.1, P6:-9.2,
    X1:14.3, X2:12.7, X3:13.9, X4:8.7,
    N1d:10.7, N1a:10.7, N2d:7.8,  N2a:7.8,  N3d:3.8,  N3a:3.8,
    N4d:4.3,  N4a:4.3,  N5d:2.2,  N5a:2.2,  N6d:1.0,  N6a:1.0,
    P1d:0.2,  P1a:0.2,  P2d:-1.9, P2a:-1.9, P3d:-3.5, P3a:-3.5,
    P4d:-5.1, P4a:-5.1, P5d:-7.0, P5a:-7.0, P6d:-7.4, P6a:-7.4,
    Q1:-10.9, Q2:-15.1, Q3:-17.4, Q4:-18.8, Q5:-23.0, D:-26.8,

    SC1:14.2, SC2:9.9,  SC3:9.2,  SC4:8.4,  SC5:6.3,  SC6:5.3,
    SN1:3.6,  SN2:2.1,  SN3:-1.8, SN4:-0.9, SN5:-3.6, SN6:-4.2,
    SP1:-5.2, SP2:-6.9, SP3:-7.7, SP4:-9.8, SP5:-11.8, SP6:-12.0,
    SX1:9.4,  SX2:7.2,  SX3:8.0,  SX4:4.3,
    SN1d:6.0,  SN1a:6.0,  SN2d:3.8,  SN2a:3.8,  SN3d:0.2,  SN3a:0.2,
    SN4d:1.1,  SN4a:1.1,  SN5d:-1.0, SN5a:-1.0, SN6d:-2.5, SN6a:-2.5,
    SP1d:-3.7, SP1a:-3.7, SP2d:-5.4, SP2a:-5.4, SP3d:-6.1, SP3a:-6.1,
    SP4d:-7.8, SP4a:-7.8, SP5d:-9.5, SP5a:-9.5, SP6d:-9.6, SP6a:-9.6,
    SQ1:-10.6, SQ2:-14.3, SQ3:-18.0, SQ4:-18.2, SQ5:-18.2, SD:-36.4,

    TC1:12.0, TC2:7.8,  TC3:6.7,  TC4:6.4,  TC5:4.5,  TC6:3.6,
    TN1:2.3,  TN2:0.3,  TN3:-3.1, TN4:-2.9, TN5:-4.9, TN6:-6.1,
    TP1:-7.2, TP2:-8.8, TP3:-9.8, TP4:-12.1, TP5:-15.2, TP6:-14.8,
    TX1:7.6,  TX2:5.2,  TX3:5.4,  TX4:2.7,
    TN1d:3.9,  TN1a:3.9,  TN2d:2.3,  TN2a:2.3,  TN3d:-1.4, TN3a:-1.4,
    TN4d:-1.2, TN4a:-1.2, TN5d:-2.8, TN5a:-2.8, TN6d:-4.1, TN6a:-4.1,
    TP1d:-5.0, TP1a:-5.0, TP2d:-6.8, TP2a:-6.8, TP3d:-7.8, TP3a:-7.8,
    TP4d:-9.5, TP4a:-9.5, TP5d:-13.2, TP5a:-13.2, TP6d:-12.7, TP6a:-12.7,
    TQ1:-14.2, TQ2:-14.5, TQ3:-18.7, TQ4:-16.3, TQ5:-17.0, TD:-36.8,
};

// Build a canonicalized version of FRAG_DELTA_F using RDKit to normalize SMILES.
// Called once per predict session. Keys in the original table may be in aromatic
// notation (e.g. "cc") or non-canonical form; this ensures the lookup always
// succeeds when comparing against RDKit-canonicalized fragment SMILES.
function buildCanonTable(RDKit) {
    const result = {};
    for (const [smi, df] of Object.entries(FRAG_DELTA_F)) {
        try {
            const mol = RDKit.get_mol(smi);
            if (!mol) {
                // Open-chain aromatic SMILES like "cc", "ccc" are invalid for
                // RDKit (no ring to close) but are legitimate AutoMartini table
                // keys. Keep them as-is for direct lookup by aromatic fragments.
                if (!(smi in result)) result[smi] = df;
                continue;
            }
            const canon = mol.get_smiles();
            mol.delete();
            if (canon && !(canon in result)) result[canon] = df;
        } catch (_) {
            // RDKit threw on an open-chain aromatic key — keep it verbatim so the
            // aromatic-notation fragment lookup can still match it.
            if (!(smi in result)) result[smi] = df;
        }
    }
    return result;
}

function _closestType(deltaF, candidates) {
    let best = candidates[0];
    let bestErr = Infinity;
    for (const t of candidates) {
        if (DELTA_F[t] === undefined) continue;
        const err = Math.abs(DELTA_F[t] - deltaF);
        if (err < bestErr) { bestErr = err; best = t; }
    }
    return best;
}

// Port of AutoMartini M3 determine_bead_type().
// deltaF: free energy of transfer (kJ/mol) — from table or Crippen conversion.
// heavyCount: non-H atoms in bead; inRing: any bead atom is aromatic.
function determineBeadType({ deltaF, charge, hDonors, hAcceptors, hasHalogen, inRing, heavyCount }) {
    const sz = heavyCount <= 2 ? 'T' : heavyCount === 3 ? 'S' : '';
    const p = sz; // size prefix

    let result;

    if (charge !== 0) {
        result = _closestType(deltaF,
            [`${p}Q1`,`${p}Q2`,`${p}Q3`,`${p}Q4`,`${p}Q5`,`${p}D`]);
    } else if (hAcceptors > 0 && hDonors === 0) {
        // Pure acceptor → 'a' label.
        result = _closestType(deltaF, p === 'T'
            ? ['TN1a','TN2a','TN3a','TN4a','TN5a','TN6a','TP1a','TP2a','TP3a','TP4a','TP5a','TP6a']
            : p === 'S'
            ? ['SN1a','SN2a','SN3a','SN4a','SN5a','SN6a','SP1a','SP2a','SP3a','SP4a','SP5a','SP6a']
            : ['N1a','N2a','N3a','N4a','N5a','N6a','P1a','P2a','P3a','P4a','P5a','P6a']);
    } else if (hDonors > 0 && hAcceptors === 0) {
        // Pure donor → 'd' label.
        result = _closestType(deltaF, p === 'T'
            ? ['TN1d','TN2d','TN3d','TN4d','TN5d','TN6d','TP1d','TP2d','TP3d','TP4d','TP5d','TP6d']
            : p === 'S'
            ? ['SN1d','SN2d','SN3d','SN4d','SN5d','SN6d','SP1d','SP2d','SP3d','SP4d','SP5d','SP6d']
            : ['N1d','N2d','N3d','N4d','N5d','N6d','P1d','P2d','P3d','P4d','P5d','P6d']);
    } else {
        // Either no H-bonding, OR both donor AND acceptor (a balanced group such
        // as an amide) → plain N/P/C series with no a/d label, matching Martini 3
        // where the a/d suffixes denote one-sided H-bonders only.
        let cands;
        if (p === 'T') {
            cands = ['TP6','TP5','TP4','TP3','TP2','TP1','TC6','TC5','TC4','TC3','TC2','TC1','TN6','TN5','TN4','TN3','TN2','TN1'];
            if (!inRing) cands = cands.filter(t => t !== 'TC5');
        } else if (p === 'S') {
            cands = ['SP6','SP5','SP4','SP3','SP2','SP1','SC6','SC5','SC4','SC3','SC2','SC1','SN6','SN5','SN4','SN3','SN2','SN1'];
        } else {
            cands = ['P6','P5','P4','P3','P2','P1','C6','C5','C4','C3','C2','C1','N6','N5','N4','N3','N2','N1'];
        }
        result = _closestType(deltaF, cands);
    }

    // Halogen overrides everything (same as AutoMartini)
    if (hasHalogen) {
        result = _closestType(deltaF, p === 'T' ? ['TX4','TX3','TX2','TX1']
                                    : p === 'S' ? ['SX4','SX3','SX2','SX1']
                                    :              ['X4','X3','X2','X1']);
    }

    return result;
}

/* ===========================================================================
   Shrake-Rupley SASA
   ===========================================================================
   Standard Bondi vdW radii (Angstrom) for elements found in organic molecules.
   Used for the AA surface numerical computation. */
// Radii in Angstrom (converted from nm). Source: the same table used by the
// reference GROMACS pipeline this tool is designed to complement.
const VDW_RADII = {
    H: 1.09, C: 1.75, N: 1.61, O: 1.56, F: 1.44,
    P: 1.80, S: 1.79, CL: 1.74, BR: 1.85, I: 2.00,
};
const DEFAULT_VDW_RADIUS = 1.75;

// Fibonacci golden-ratio lattice — generates n near-uniformly distributed
// points on the unit sphere.
function fibonacciSpherePoints(n) {
    const pts = new Array(n);
    const phi = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < n; i++) {
        const y = 1 - (i / (n - 1)) * 2;
        const r = Math.sqrt(1 - y * y);
        const theta = phi * i;
        pts[i] = [r * Math.cos(theta), y, r * Math.sin(theta)];
    }
    return pts;
}

// Shrake-Rupley numerical SASA in Å².
// particles: array of [x, y, z, vdwRadius] in Angstrom.
// nPoints: number of test points per sphere; 960 gives ~1% accuracy for typical
// drug-like molecules and runs in milliseconds.
function shrakeRupley(particles, probeRadius, nPoints = 4800) {
    const n = particles.length;
    if (n === 0) return 0;
    const unitPts = fibonacciSpherePoints(nPoints);
    let totalSASA = 0;

    for (let i = 0; i < n; i++) {
        const [xi, yi, zi, ri] = particles[i];
        const shellR = ri + probeRadius;

        // Pre-filter: j can only bury points on i's shell if centers are closer
        // than shellR_i + shellR_j. Building this list once per i is much faster
        // than checking all N for each of the nPoints test points.
        const neighbors = [];
        for (let j = 0; j < n; j++) {
            if (j === i) continue;
            const [xj, yj, zj, rj] = particles[j];
            const cutoff = shellR + rj + probeRadius;
            const dx = xi - xj, dy = yi - yj, dz = zi - zj;
            if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) {
                neighbors.push(j);
            }
        }

        let exposed = 0;
        for (const [ux, uy, uz] of unitPts) {
            const px = xi + shellR * ux;
            const py = yi + shellR * uy;
            const pz = zi + shellR * uz;
            let buried = false;
            for (const j of neighbors) {
                const [xj, yj, zj, rj] = particles[j];
                const cutoff = rj + probeRadius;
                const dx = px - xj, dy = py - yj, dz = pz - zj;
                if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) {
                    buried = true;
                    break;
                }
            }
            if (!buried) exposed++;
        }
        totalSASA += (exposed / nPoints) * 4 * Math.PI * shellR * shellR;
    }
    return totalSASA;
}

// Extract per-atom positions and vdW radii from an NGL Structure, then compute
// SASA with Shrake-Rupley. Hydrogen atoms contribute negligibly; they are
// included to match standard SASA conventions.
function aaSASA(structure, probeRadius) {
    const particles = [];
    structure.eachAtom((atom) => {
        const el = (atom.element || "").toUpperCase();
        const r = VDW_RADII[el] ?? DEFAULT_VDW_RADIUS;
        particles.push([atom.x, atom.y, atom.z, r]);
    });
    return shrakeRupley(particles, probeRadius);
}

// Compute CG SASA from the bead collection using Martini bead radii.
// Virtual (U) beads and empty beads are excluded.
function cgSASA(collection, probeRadius) {
    const particles = [];
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) continue;
        const r = beadRadius(bead);
        if (r <= 0) continue;
        const c = bead.center;
        particles.push([c.x, c.y, c.z, r]);
    }
    return shrakeRupley(particles, probeRadius);
}

/* ===========================================================================
   Fragment SMILES builder
   ===========================================================================
   Generates a SMILES string for the subgraph of atoms belonging to one bead,
   capping open valences (bonds to atoms outside the bead) with implicit H.
   Handles rings via DFS back-edge detection. */

const _ATOM_VALENCE = {C:4,N:3,O:2,S:2,P:3,F:1,CL:1,BR:1,I:1,H:1,SE:2,SI:4,B:3};
const _ORGANIC = new Set(['B','C','N','O','P','S','F','CL','BR','I']);

// Reference bond lengths (Å) as [single, double, triple]; null = not applicable.
// Used to infer bond order from 3D geometry when the input file carries no
// bond-order metadata (the common case for PDB/GRO). Keyed by element pair
// sorted alphabetically.
const _BOND_REF = {
    'C-C': [1.54, 1.34, 1.20],
    'C-N': [1.47, 1.29, 1.16],
    'C-O': [1.43, 1.23, null],
    'C-S': [1.82, 1.60, null],
    'C-P': [1.84, 1.66, null],
    'N-N': [1.45, 1.25, 1.10],
    'N-O': [1.41, 1.21, null],
    'N-S': [1.71, 1.54, null],
    'N-P': [1.65, 1.55, null],
    'O-O': [1.48, 1.21, null],
    'O-P': [1.63, 1.50, null],
    'O-S': [1.57, 1.44, null],
    'P-S': [2.12, 1.95, null],
};

// Classify a bond's order from its measured length. Returns {order, single}
// where `single` is the single-bond reference (used to score "how shortened").
function _geometryBondOrder(elA, elB, dist) {
    const ref = _BOND_REF[[elA, elB].sort().join('-')];
    if (!ref || !dist) return { order: 1, single: ref ? ref[0] : Infinity };
    let best = 1, bestErr = Infinity;
    for (let o = 1; o <= 3; o++) {
        const r = ref[o - 1];
        if (r == null) continue;
        const err = Math.abs(dist - r);
        if (err < bestErr) { bestErr = err; best = o; }
    }
    return { order: best, single: ref[0] };
}

// Perceive whole-molecule chemistry from geometry (and any NGL bond orders that
// happen to be present). Returns:
//   aromaticAtoms : Set of global atom indices in an aromatic ring
//   bondOrders    : Map "a-b" (a<b) -> order, for NON-aromatic multiple bonds
//
// Bond orders come from comparing each heavy-heavy bond length to reference
// single/double/triple lengths. Aromaticity = ring topology (prune + cycles)
// where every ring member is sp2 (carries a shortened/multiple bond) or a
// 2-connected heteroatom. Non-aromatic double/triple bonds are assigned greedily
// (most-shortened first, at most one multiple bond per atom) so resonance/amide
// carbons don't get an impossible double-double valence (C=O wins, C-N stays
// single). Aromatic bonds are left to the per-fragment kekuliser.
function perceiveChemistry(structure) {
    const aromaticAtoms = new Set();
    const bondOrders = new Map();
    if (!structure || typeof structure.eachAtom !== 'function') {
        return { aromaticAtoms, bondOrders };
    }

    // Element + position per atom. NGL reuses one AtomProxy during iteration, so
    // read values into plain objects immediately rather than storing the proxy.
    const element = new Map();
    const pos = new Map();
    structure.eachAtom(a => {
        element.set(a.index, (a.element || 'C').toUpperCase());
        const v = typeof a.positionToVector3 === 'function'
            ? a.positionToVector3() : { x: a.x, y: a.y, z: a.z };
        pos.set(a.index, { x: v.x, y: v.y, z: v.z });
    });

    // Heavy-atom graph + per-bond geometry classification (each bond once).
    const adj = new Map();
    const bonds = [];          // { a, b, order, single, aromaticNGL }
    const seen = new Set();
    const sp2 = new Set();
    const link = (u, v) => {
        if (!adj.has(u)) adj.set(u, new Set());
        adj.get(u).add(v);
    };
    structure.eachAtom(atom => {
        const idx = atom.index;
        if (element.get(idx) === 'H') return;
        if (atom.aromatic) sp2.add(idx);  // trust NGL's flag when present
        if (typeof atom.eachBond !== 'function') return;
        atom.eachBond(bond => {
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const other = i1 === idx ? i2 : (i2 === idx ? i1 : -1);
            if (other < 0 || element.get(other) === 'H') return;
            link(idx, other); link(other, idx);
            const a = Math.min(idx, other), b = Math.max(idx, other);
            const key = `${a}-${b}`;
            if (seen.has(key)) return;
            seen.add(key);
            const pa = pos.get(a), pb = pos.get(b);
            const dist = Math.hypot(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z);
            let { order, single } = _geometryBondOrder(element.get(a), element.get(b), dist);
            // Respect explicit NGL bond orders when the file does provide them.
            const ngl = bond.bondOrder;
            if (ngl === 2 || ngl === 3) order = Math.max(order, ngl);
            const aromaticNGL = ngl === 1.5 || (ngl >= 4 && ngl < 100);
            if (order >= 2 || aromaticNGL) { sp2.add(a); sp2.add(b); }
            bonds.push({ a, b, order, single, dist, aromaticNGL });
        });
    });

    // Prune terminal chains: atoms with degree <= 1 cannot be in a ring.
    const degree = new Map();
    for (const [k, s] of adj) degree.set(k, s.size);
    const removed = new Set();
    let changed = true;
    while (changed) {
        changed = false;
        for (const [k, d] of degree) {
            if (removed.has(k) || d > 1) continue;
            removed.add(k);
            changed = true;
            for (const nb of adj.get(k)) {
                if (!removed.has(nb)) degree.set(nb, degree.get(nb) - 1);
            }
        }
    }
    const ringAdj = new Map();
    for (const k of adj.keys()) {
        if (removed.has(k)) continue;
        ringAdj.set(k, [...adj.get(k)].filter(n => !removed.has(n)));
    }

    // Fundamental cycles via recursive DFS (molecules are small).
    const color = new Map(); // 1 = on stack, 2 = done
    const par = new Map();
    const cycles = [];
    const dfs = (u, p) => {
        color.set(u, 1); par.set(u, p);
        for (const v of ringAdj.get(u)) {
            if (v === p) continue;
            if (color.get(v) === 1) {
                const cyc = [u];
                let x = u;
                while (x !== v) { x = par.get(x); cyc.push(x); }
                cycles.push(cyc);
            } else if (!color.has(v)) {
                dfs(v, u);
            }
        }
        color.set(u, 2);
    };
    for (const start of ringAdj.keys()) {
        if (!color.has(start)) dfs(start, -1);
    }

    const isHetero = idx => ['N', 'O', 'S'].includes(element.get(idx));
    for (const cyc of cycles) {
        if (cyc.length < 4 || cyc.length > 7) continue;
        if (cyc.every(idx => sp2.has(idx) || isHetero(idx))) {
            for (const idx of cyc) aromaticAtoms.add(idx);
        }
    }

    // Assign non-aromatic multiple bonds greedily: most-shortened first, at most
    // one per atom. Aromatic bonds are excluded — kekulisation handles those.
    // The one-per-atom cap stops resonance/amide carbons getting an impossible
    // double-double valence: e.g. for an amide the C=O (more shortened) wins and
    // the C–N is left single.
    const ranked = bonds
        .filter(bd => bd.order >= 2 && !(aromaticAtoms.has(bd.a) && aromaticAtoms.has(bd.b)))
        .sort((x, y) => (y.single - y.dist) - (x.single - x.dist));
    const usedMul = new Set();
    for (const bd of ranked) {
        if (usedMul.has(bd.a) || usedMul.has(bd.b)) continue;
        bondOrders.set(`${bd.a}-${bd.b}`, bd.order);
        usedMul.add(bd.a); usedMul.add(bd.b);
    }
    return { aromaticAtoms, bondOrders };
}

// aromaticNotation=true: write aromatic atoms as lowercase (c, n, o…), skip
// kekulization. Used to generate AutoMartini-style lookup keys ("cc", "ccc"…)
// that match open-chain aromatic entries in the table. Never feed this output
// to RDKit — it rejects open-chain aromatic SMILES.
//
// aromaticAtoms: optional Set of global atom indices known to be aromatic (from
// perceiveChemistry). When provided it overrides the local flag/bond-order
// heuristic, which is unreliable for Kekulé-form input.
// bondOrders: optional Map "a-b"(a<b) -> order for non-aromatic multiple bonds
// (from perceiveChemistry's geometry inference). Overrides NGL bond orders so
// carbonyls/imines/alkenes render correctly even from PDB/GRO without bond data.
function fragmentToSmiles(beadAtoms,
        { aromaticNotation = false, aromaticAtoms = null, bondOrders = null,
          startIndex = null } = {}) {
    const heavy = beadAtoms.filter(a => (a.element || 'C').toUpperCase() !== 'H');
    if (heavy.length === 0) return null;

    const atomSet = new Set(heavy.map(a => a.index));

    // Determine which bead atoms are aromatic. Prefer the caller-supplied set
    // (whole-molecule ring perception); otherwise fall back to the local
    // flag/bond-order heuristic (works only when NGL encodes aromaticity).
    const aromaticSet = new Set();
    if (aromaticAtoms) {
        for (const atom of heavy) {
            if (aromaticAtoms.has(atom.index)) aromaticSet.add(atom.index);
        }
    } else {
        for (const atom of heavy) {
            if (atom.aromatic) { aromaticSet.add(atom.index); continue; }
            if (typeof atom.eachBond !== 'function') continue;
            let found = false;
            atom.eachBond(b => {
                const bo = b.bondOrder;
                if (bo === 1.5 || (bo >= 4 && bo < 100)) found = true;
            });
            if (found) aromaticSet.add(atom.index);
        }
    }

    // Collect each internal bond once (undirected), flagging aromatic bonds.
    // NGL may encode aromaticity as bond order 4 or 1.5, or only via the per-atom
    // aromatic flag, so check both.
    const edges = new Map(); // "min-max" -> { a, b, order, aromatic }
    for (const atom of heavy) {
        if (typeof atom.eachBond !== 'function') continue;
        atom.eachBond(bond => {
            // Identify the other endpoint robustly: NGL may order atomIndex1/2
            // with the lower index first (not necessarily the current atom first).
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const otherIdx = i1 === atom.index ? i2 : (i2 === atom.index ? i1 : -1);
            if (otherIdx < 0 || !atomSet.has(otherIdx)) return;
            const a = Math.min(atom.index, otherIdx);
            const b = Math.max(atom.index, otherIdx);
            const key = `${a}-${b}`;
            if (edges.has(key)) return;
            const raw = bond.bondOrder || 1;
            const aromatic = raw >= 4 || raw === 1.5
                || (aromaticSet.has(a) && aromaticSet.has(b));
            // Bond order: aromatic bonds get kekulised below; otherwise prefer
            // the geometry-inferred order, falling back to NGL's bond order.
            let order;
            if (aromatic) order = 1;
            else if (bondOrders && bondOrders.has(key)) order = bondOrders.get(key);
            else order = Math.round(raw);
            edges.set(key, { a, b, order, aromatic });
        });
    }

    // Kekulise aromatic bonds (Kekulé mode only): a greedy maximum matching
    // turns alternating aromatic bonds into double bonds so sp2 atoms keep
    // their unsaturation. Skipped in aromaticNotation mode where aromatic atoms
    // are written lowercase and bonds are implicit (order stays 1).
    if (!aromaticNotation) {
        const doubled = new Set();
        for (const edge of edges.values()) {
            if (edge.aromatic && !doubled.has(edge.a) && !doubled.has(edge.b)) {
                edge.order = 2;
                doubled.add(edge.a);
                doubled.add(edge.b);
            }
        }
    }

    // Per-atom data: internal bonds and implicit H count. Any valence not used by
    // an internal bond is capped with hydrogen (external heavy bonds and original
    // explicit hydrogens alike).
    const data = new Map();
    for (const atom of heavy) {
        data.set(atom.index, {
            el: (atom.element || 'C').toUpperCase(),
            charge: atom.formalCharge ?? 0,
            hCount: 0,
            internalBonds: [],
        });
    }
    const intSum = new Map();
    for (const edge of edges.values()) {
        data.get(edge.a).internalBonds.push({ toIdx: edge.b, order: edge.order });
        data.get(edge.b).internalBonds.push({ toIdx: edge.a, order: edge.order });
        intSum.set(edge.a, (intSum.get(edge.a) || 0) + edge.order);
        intSum.set(edge.b, (intSum.get(edge.b) || 0) + edge.order);
    }
    for (const atom of heavy) {
        const d = data.get(atom.index);
        const baseVal = _ATOM_VALENCE[d.el] ?? 4;
        d.hCount = Math.max(0, baseVal + d.charge - (intSum.get(atom.index) || 0));
    }

    // DFS to find ring-closure back-edges and assign ring-closure digits
    const visited = new Set();
    const inStack = new Set();
    const seenEdges = new Set();
    const ringClosures = new Map(); // atomIdx -> [{digit, writeBond}]
    let digitCounter = 1;

    function findBackEdges(idx, parentIdx) {
        visited.add(idx); inStack.add(idx);
        for (const { toIdx, order } of (data.get(idx)?.internalBonds ?? [])) {
            const edgeKey = `${Math.min(idx, toIdx)}-${Math.max(idx, toIdx)}`;
            if (toIdx === parentIdx || seenEdges.has(edgeKey)) continue;
            seenEdges.add(edgeKey);
            if (inStack.has(toIdx)) {
                // Back-edge: toIdx is the opener (visited first), idx is the closer
                const digit = digitCounter++;
                if (!ringClosures.has(toIdx)) ringClosures.set(toIdx, []);
                if (!ringClosures.has(idx))   ringClosures.set(idx,   []);
                ringClosures.get(toIdx).push({ digit, writeBond: false });
                ringClosures.get(idx).push({ digit, order, writeBond: true });
            } else if (!visited.has(toIdx)) {
                findBackEdges(toIdx, idx);
            }
        }
        inStack.delete(idx);
    }

    // DFS start atom: caller may pin it (to enumerate alternative serialisations
    // for table matching); default to the first heavy atom.
    const startIdx = (startIndex != null && atomSet.has(startIndex))
        ? startIndex : heavy[0].index;
    findBackEdges(startIdx, -1);
    visited.clear();

    function bondChar(order) {
        return order === 2 ? '=' : order === 3 ? '#' : '';
    }

    function atomToken(idx) {
        const d = data.get(idx);
        // In aromatic mode, write aromatic atoms as lowercase (c, n, o…).
        // Bonds between lowercase atoms are implicitly aromatic so no bond
        // character is needed — bondChar(1) already returns ''.
        if (aromaticNotation && aromaticSet.has(idx)) {
            return d.el.charAt(0).toLowerCase() + d.el.slice(1).toLowerCase();
        }
        const sym = d.el === 'CL' ? 'Cl' : d.el === 'BR' ? 'Br'
                  : d.el.charAt(0) + d.el.slice(1).toLowerCase();

        const needsBracket = d.hCount > 0 || d.charge !== 0 || !_ORGANIC.has(d.el);
        if (!needsBracket) return sym;

        let inner = sym;
        if (d.hCount === 1) inner += 'H';
        else if (d.hCount > 1) inner += `H${d.hCount}`;
        if (d.charge > 0) inner += d.charge === 1 ? '+' : `+${d.charge}`;
        else if (d.charge < 0) inner += d.charge === -1 ? '-' : `${d.charge}`;
        return `[${inner}]`;
    }

    function closureSuffix(idx) {
        return (ringClosures.get(idx) ?? []).map(({ digit, order, writeBond }) => {
            const b = writeBond ? bondChar(order) : '';
            return b + (digit >= 10 ? `%${digit}` : `${digit}`);
        }).join('');
    }

    function dfs(idx) {
        visited.add(idx);
        let smi = atomToken(idx) + closureSuffix(idx);

        const children = (data.get(idx)?.internalBonds ?? [])
            .filter(({ toIdx }) => !visited.has(toIdx));

        if (children.length === 0) return smi;
        // Branches must immediately follow the parent atom. Write every child
        // except the last as a parenthesised branch, then the last child inline
        // as the main-chain continuation. (Appending branches AFTER an inline
        // child would wrongly attach them to the tail of that child's subtree —
        // e.g. an amide carbon C(N)(O) became C-N-O.) Children visited during a
        // sibling's recursion (rings) are skipped here; their ring-closure digit
        // is emitted via closureSuffix instead.
        for (let i = 0; i < children.length - 1; i++) {
            const { toIdx, order } = children[i];
            if (visited.has(toIdx)) continue;
            smi += `(${bondChar(order)}${dfs(toIdx)})`;
        }
        const last = children[children.length - 1];
        if (!visited.has(last.toIdx)) {
            smi += bondChar(last.order) + dfs(last.toIdx);
        }
        return smi;
    }

    return dfs(startIdx);
}

/* ===========================================================================
   H-bond donor detection (from the intact structure)
   ===========================================================================
   Donor status must be read from the real molecule, NOT the capped fragment
   SMILES: capping a heteroatom's broken bonds with hydrogen turns non-donors
   (tertiary amines, ethers) into false donors. Acceptor status survives capping
   (lone pairs are preserved), so that stays with RDKit's NumHBA. */

const _HETERO_VALENCE = { N: 3, O: 2, S: 2 };

// Number of hydrogens actually attached to a heteroatom in the real structure.
// hasExplicitH: whether the loaded molecule carries explicit hydrogen atoms.
function realHydrogenCount(atom, hasExplicitH) {
    const el = (atom.element || '').toUpperCase();
    const baseVal = _HETERO_VALENCE[el];
    if (baseVal === undefined) return 0;

    const charge = atom.formalCharge ?? 0;
    const structure = atom.structure;
    let explicitH = 0;
    let heavySum = 0;
    let aromaticRingBonds = 0;

    if (typeof atom.eachBond === 'function') {
        atom.eachBond((bond) => {
            const otherIdx = bond.atomIndex1 === atom.index
                ? bond.atomIndex2 : bond.atomIndex1;
            const other = structure.getAtomProxy(otherIdx);
            if ((other.element || '').toUpperCase() === 'H') { explicitH += 1; return; }
            let order = bond.bondOrder || 1;
            // Treat as an aromatic ring bond if both atoms are flagged aromatic,
            // OR the bond order encodes aromaticity (NGL may use 4 or 1.5).
            const aromatic = (!!atom.aromatic && !!other.aromatic)
                || order >= 4 || order === 1.5;
            if (aromatic) { aromaticRingBonds += 1; order = 1; }
            heavySum += order;
        });
    }

    // If the molecule carries explicit hydrogens, the count of H neighbours is
    // definitive — even when it is zero (e.g. a pyridine-type aromatic N).
    if (hasExplicitH) return explicitH;

    // No explicit H: estimate implicit H from valence. An aromatic ring atom
    // consumes one extra valence unit beyond its (normalised) single ring bonds,
    // so pyridine N -> 0 H and benzene C -> 1 H. (A pyrrole-type N-H cannot be
    // told apart from a pyridine N without explicit hydrogens.)
    let used = heavySum;
    if (aromaticRingBonds >= 2) used += 1;
    return Math.max(0, baseVal + charge - used);
}

// Count heteroatoms in a bead that bear at least one hydrogen (H-bond donors).
function beadDonorCount(bead, hasExplicitH) {
    let donors = 0;
    for (const atom of bead.atoms) {
        const el = (atom.element || '').toUpperCase();
        if ((el === 'N' || el === 'O' || el === 'S')
            && realHydrogenCount(atom, hasExplicitH) > 0) {
            donors += 1;
        }
    }
    return donors;
}

// True if the loaded NGL structure contains any explicit hydrogen atoms.
function structureHasHydrogens(structure) {
    let found = false;
    if (structure && typeof structure.eachAtom === 'function') {
        structure.eachAtom((atom) => {
            if ((atom.element || '').toUpperCase() === 'H') found = true;
        });
    }
    return found;
}

/* ===========================================================================
   RDKit.js lazy loader
   =========================================================================== */

let _rdkitPromise = null;

function loadRDKit() {
    if (_rdkitPromise) return _rdkitPromise;
    _rdkitPromise = new Promise((resolve, reject) => {
        if (window.RDKit) { resolve(window.RDKit); return; }
        const script = document.createElement('script');
        script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
        script.onload = () => {
            window.initRDKitModule()
                .then(rdkit => { window.RDKit = rdkit; resolve(rdkit); })
                .catch(reject);
        };
        script.onerror = () => reject(new Error('Failed to load RDKit.js'));
        document.head.appendChild(script);
    });
    return _rdkitPromise;
}

function loadMolecule(event, stage) {
    // Clear the stage if needed
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();
    // Setup the model
    let collection = new BeadCollection();
    // Setup the interface
    let vizu = new Visualization(collection, stage);
    // Load the molecule
    let input = event.target.files[0];

    const namePromise = readOriginalAtomNames(input);
    const componentPromise = stage.loadFile(input);

    Promise.all([namePromise, componentPromise])
        .then(([names, component]) => {
            // Ensure original atom names are applied before any UI/render
            collection.setOriginalAtomNames(names || []);
            component.addRepresentation("ball+stick", bondAwareRepresentationParams());
            component.autoView();
            vizu.attachAALabels(component);
            vizu.attachRepresentation(component);
            vizu.updateSelection();
        })
        .catch((err) => {
            console.error("Error loading molecule or reading original atom names:", err);
        });
    // Bing the new bead buttons.
    let buttons = document.getElementsByClassName("new-bead");
    for (const button of buttons) {
        button.onclick = (event) => vizu.onNewBead(event);
        button.disabled = false;
    }
	// Bind our own selection behaviour.
    // We need to use the "arrow" function so that `this` is defined and refer
    // to the right object in the `onClick` method. See
    // <https://stackoverflow.com/questions/20279484/how-to-access-the-correct-this-inside-a-callback>.
    stage.signals.clicked.add((pickingProxy) => vizu.onClick(pickingProxy));
}

function main() {
    // Capture the wheel events within the viewer so the page does not scroll when we zoom in or out.
    // <https://github.com/nglviewer/ngl/issues/878#issuecomment-913504711>
    const stageContainer = document.getElementById('viewport');
    function maybeScroll(event) {
        if (stageContainer.contains(event.target)) {     // If wheel event occurred within the viewer
            event.preventDefault();                      // prevent the default (scrolling the page)
        }
    }
    window.addEventListener('wheel', maybeScroll, {passive: false});

    // Create NGL Stage object
    let stage = new NGL.Stage( "viewport" );

    // Handle window resizing
    window.addEventListener( "resize", function( event ){
        stage.handleResize();
    }, false );

	let mol_select = document.getElementById("mol-select");
	mol_select.onchange = (event) => loadMolecule(event, stage);
	
	// Remove preset action on atom pick.
	// As of NGL v2.0.0-dev.11, the left click atom pick is bind to the
	// centering of the view on the selected atom. In previous versions, this
	// behavior was linked on shift-click, instead.
	stage.mouseControls.remove("clickPick-left");

    let buttons = document.getElementsByClassName("new-bead");
    for (const button of buttons) {
        button.disabled = true;
    }
}

window.onload = main;