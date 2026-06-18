import { buildCanonTable, determineBeadType } from './prediction.js';
import { perceiveChemistry, fragmentToSmiles, beadDonorCount, structureHasHydrogens } from './chemistry.js';
import { PROBE_RADIUS, aaSASA, cgSASA, beadsToPDB } from './sasa.js';
import { generateNDX, generateMap, generateGRO, generatePythonAssignments,
         download, copyTextToClipboard, bondAwareRepresentationParams,
         parseShakerMapping } from './fileformats.js';
import { loadRDKit } from './rdkit.js';

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

function findParentWithClass(element, className) {
    let node = element;
    while (node) {
        if (node.classList.contains(className)) return node;
        node = node.parentElement;
    }
    return null;
}

export class Visualization {
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
                setTimeout(() => { button.textContent = originalText; }, 1200);
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
    }

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

    onToggleCG(event) {
        this.showCG = event.target.checked;
        this.drawCG();
    }

    onToggleAASurface(event) {
        this.showAASurface = event.target.checked;
        if (this.aaSurface) this.aaSurface.setVisibility(this.showAASurface);
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

            const canonTable = buildCanonTable(RDKit);

            btn.textContent = 'Predicting…';
            const hasExplicitH = structureHasHydrogens(
                this.component && this.component.structure);
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

    onToggleAALabels(event) {
        this.aa_labels.setVisibility(event.target.checked);
    }

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

	onNewBead(event) {
	    this.collection.newBead();
	    this.updateSelection();
	}

    onBeadSelected(event) {
        const tag = event.target.tagName;
        if (tag === "INPUT" || tag === "BUTTON" || tag === "FORM" || tag === "LABEL") return;

        let realTarget = findParentWithClass(event.target, "bead-view");
        let nodes = document.getElementById("bead-list").childNodes;
        let index = 0;
        for (const child of nodes) {
            if (child === realTarget) {
                if (this.collection.beads[index] === this.currentBead) {
                    this.collection.deselectBead();
                } else {
                    this.collection.selectBead(index);
                }
                break;
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
            if (child === realTarget) { selected = index; break; }
            index += 1;
        }
        if (selected >= 0) {
            this.collection.removeBead(selected);
            if (this.collection.beads.length === 0) this.collection.newBead();
            if (realTarget.classList.contains('selected-bead')) this.collection.selectBead(0);
        }
        this.updateSelection();
    }

    onNameChange(event) {
        let realTarget = findParentWithClass(event.target, "bead-view");
        let nodes = document.getElementById("bead-list").childNodes;
        let index = 0;
        for (const child of nodes) {
            if (child === realTarget) this.collection.beads[index].name = event.target.value;
            index += 1;
        }
        this.updateName();
        this.checkDuplicateNames();
    }

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
        addLabeledField("Charge", chargeNode);

        // DELETE BUTTON
        let removeNode = document.createElement("button");
        removeNode.textContent = "Delete";
        removeNode.classList.add("delete-bead");
        removeNode.onclick = (event) => { event.stopPropagation(); this.onBeadRemove(event); };

        headerRow.appendChild(fieldsNode);
        headerRow.appendChild(removeNode);
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

    createBeadList() {
        const counts = new Map();
        for (const bead of this.collection.beads)
            counts.set(bead.name, (counts.get(bead.name) || 0) + 1);
        const dupes = new Set([...counts.entries()].filter(([, n]) => n > 1).map(([k]) => k));
        for (const bead of this.collection.beads) this.createBeadListItem(bead, dupes.has(bead.name));
    }

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

    clearBeadList() {
        let list = document.getElementById('bead-list');
        while (list.lastChild) list.removeChild(list.lastChild);
    }

    updateNDX() { document.getElementById('ndx-output').textContent = generateNDX(this.collection); }
    updateMap() { document.getElementById('map-output').textContent = generateMap(this.collection); }
    updateGRO() { document.getElementById('gro-output').textContent = generateGRO(this.collection); }
    updatePY() {
        const counts = new Map();
        for (const bead of this.collection.beads)
            counts.set(bead.name, (counts.get(bead.name) || 0) + 1);
        const hasDupes = [...counts.values()].some(n => n > 1);
        document.getElementById('py-warning').hidden = !hasDupes;
        document.getElementById('py-output').textContent = generatePythonAssignments(this.collection);
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

    drawCG() {
        let selectedColor = [0.25, 0.84, 0.96];
        let opacity = this.showCG ? 1 : 0.2;

        if (this.shapeComp != null) this.stage.removeComponent(this.shapeComp);

        let shape = new NGL.Shape("shape");
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
