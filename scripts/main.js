import { BeadCollection } from './bead.js';
import { Visualization } from './visualization.js';
import { readOriginalAtomNames, bondAwareRepresentationParams } from './fileformats.js';

function loadMolecule(event, stage) {
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();

    let collection = new BeadCollection();
    let vizu = new Visualization(collection, stage);
    let input = event.target.files[0];

    const namePromise = readOriginalAtomNames(input);
    const componentPromise = stage.loadFile(input);

    Promise.all([namePromise, componentPromise])
        .then(([names, component]) => {
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

    let buttons = document.getElementsByClassName("new-bead");
    for (const button of buttons) {
        button.onclick = (event) => vizu.onNewBead(event);
        button.disabled = false;
    }

    stage.signals.clicked.add((pickingProxy) => vizu.onClick(pickingProxy));
}

function main() {
    // Capture wheel events within the viewer so the page doesn't scroll when zooming.
    // https://github.com/nglviewer/ngl/issues/878#issuecomment-913504711
    const stageContainer = document.getElementById('viewport');
    window.addEventListener('wheel', (event) => {
        if (stageContainer.contains(event.target)) event.preventDefault();
    }, { passive: false });

    let stage = new NGL.Stage("viewport");

    window.addEventListener("resize", () => stage.handleResize(), false);

	let mol_select = document.getElementById("mol-select");
	mol_select.onchange = (event) => loadMolecule(event, stage);

    // Remove preset left-click centering behaviour (added in NGL v2.0.0-dev.11).
	stage.mouseControls.remove("clickPick-left");

    let buttons = document.getElementsByClassName("new-bead");
    for (const button of buttons) button.disabled = true;
}

window.onload = main;
