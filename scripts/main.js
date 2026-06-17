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

	toggleAtom(atom) {
	    if (this.isAtomIn(atom)) {
	        this.removeAtom(atom);
	    } else {
	        this.addAtom(atom);
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
        addLabeledField("Type", typeNode);

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
            item.scrollIntoView(false);
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