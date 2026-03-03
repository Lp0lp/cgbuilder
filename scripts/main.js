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
}


class Visualization {
    constructor(collection, stage) {
        this.collection = collection;
        this.representation = null;
        this.stage = stage;
        this.shapeComp = null;
        this.showCG = false;

        let toggleCG = document.getElementById('toggle-cg');
        toggleCG.onclick = (event) => this.onToggleCG(event);
        toggleCG.disabled = false;

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
    }

	get currentBead() {
	    return this.collection.currentBead;
	}

    attachRepresentation(component) {
        this.representation = component.addRepresentation(
	        "ball+stick",
	        {
	            sele: "not all",
	            radiusScale: 1.6,
	            color: "#f4b642",
	            opacity: 0.6
	        },
	    );
    }

    attachAALabels(component) {
        this.aa_labels = component.addRepresentation(
            "label",
            {
                labelType: "atomname",
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
                const name = atom.atomname;

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
            atomname = atom.atomname;
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
        output += index + "\t" + atom.atomname;
        for (const bead of atomToBeads[atom.atomname]) {
            output += "\t" + bead;
        }
        output += "\n";
    }

    return output;
}

function generatePythonAssignments(collection) {
    const beads = collection.beads;

    if (beads.length === 0) {
        return "";
    }

    let lines = [];
    let beadVarNames = [];
    let beadTypes = [];
    let beadCharges = [];

    const resname = collection.beads[0].resname;
    lines.push(`resname='${resname}'`);

    for (const bead of beads) {

        const varName = bead.name;
        beadVarNames.push(varName);

        const atomNames = bead.expandedAtoms
            ? bead.expandedAtoms().map(a => `'${a.atomname}'`).join(",")
            : bead.atoms.map(a => `'${a.atomname}'`).join(",");

        lines.push(`${varName}   = [${atomNames}]`);

        beadTypes.push(`'${bead.type || ""}'`);
        beadCharges.push(bead.charge ?? 0);
    }

    lines.push("");
    lines.push(`bead_assignments = [${beadVarNames.join(",       ")}]`);
    lines.push(`bead_types       = [ ${beadTypes.join(", ")} ]`);
    lines.push(`bead_names       = [ ${beadVarNames.map(n => `'${n}'`).join(", ")} ]`);
    lines.push(`bead_charges     = [ ${beadCharges.join(", ")} ]`);

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


function loadMolecule(event, stage) {
    // Clear the stage if needed
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();
    // Setup the model
    let collection = new BeadCollection();
    // Setup the interface
    let vizu = new Visualization(collection, stage);
    // Load the molecule
    let input = event.target.files[0]
	stage.loadFile(input).then(function (component) {
	    component.addRepresentation("ball+stick");
	    component.autoView();
	    vizu.attachAALabels(component);
	    vizu.attachRepresentation(component);
	    vizu.updateSelection();
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