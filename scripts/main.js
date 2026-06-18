import { BeadCollection } from './bead.js';
import { Visualization } from './visualization.js';
import { readOriginalAtomNames, bondAwareRepresentationParams } from './fileformats.js';
import { EXAMPLE_PDB, EXAMPLE_MAPPING } from './example.js';

let currentVizu = null;

function loadMoleculeFromFile(file, stage) {
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();

    let collection = new BeadCollection();
    let vizu = new Visualization(collection, stage);
    currentVizu = vizu;

    const namePromise = readOriginalAtomNames(file);
    const componentPromise = stage.loadFile(file);

    const ready = Promise.all([namePromise, componentPromise])
        .then(([names, component]) => {
            collection.setOriginalAtomNames(names || []);
            component.addRepresentation("ball+stick", bondAwareRepresentationParams());
            component.autoView();
            vizu.attachAALabels(component);
            vizu.attachRepresentation(component);
            vizu.checkAtomNameUniqueness(component.structure);
            vizu.countHeavyAtoms(component.structure);
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

    document.getElementById('load-mapping-btn').disabled = false;
    document.getElementById('clear-beads-btn').disabled = false;

    stage.signals.clicked.add((pickingProxy) => vizu.onClick(pickingProxy));

    return ready;
}

function loadMolecule(event, stage) {
    loadMoleculeFromFile(event.target.files[0], stage);
}

function initTabs() {
    const btns = document.querySelectorAll('.tab-btn');
    const panels = document.querySelectorAll('.tab-panel');
    btns.forEach(btn => {
        btn.onclick = () => {
            btns.forEach(b => b.classList.remove('active'));
            panels.forEach(p => p.classList.remove('active'));
            btn.classList.add('active');
            document.getElementById('tab-' + btn.dataset.tab).classList.add('active');
        };
    });
}

function initTheme(stage) {
    const themeBtn = document.getElementById('theme-toggle');
    const bgToggle = document.getElementById('toggle-bg');

    function syncBg(isDark) {
        bgToggle.checked = !isDark;
        stage.setParameters({ backgroundColor: isDark ? 'black' : 'white' });
    }

    // Initialise BG to match the starting theme (set by the inline <head> script)
    syncBg(document.documentElement.getAttribute('data-theme') === 'dark');

    themeBtn.onclick = () => {
        const nextDark = document.documentElement.getAttribute('data-theme') !== 'dark';
        if (nextDark) {
            document.documentElement.setAttribute('data-theme', 'dark');
            localStorage.setItem('cgbuilder-theme', 'dark');
        } else {
            document.documentElement.removeAttribute('data-theme');
            localStorage.setItem('cgbuilder-theme', 'light');
        }
        syncBg(nextDark);
    };

    bgToggle.onchange = () => {
        stage.setParameters({ backgroundColor: bgToggle.checked ? 'white' : 'black' });
    };
}

function initNavbar(stage) {
    const tabs = document.querySelectorAll('.navbar-tab');
    const pages = document.querySelectorAll('.page');
    tabs.forEach(tab => {
        tab.onclick = () => {
            tabs.forEach(t => t.classList.remove('active'));
            pages.forEach(p => p.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById('page-' + tab.dataset.page).classList.add('active');
            if (tab.dataset.page === 'app') stage.handleResize();
        };
    });
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

    document.getElementById('load-example-btn').onclick = () => {
        const file = new File([EXAMPLE_PDB], 'example_L1.pdb', { type: 'text/plain' });
        loadMoleculeFromFile(file, stage).then(() => {
            if (currentVizu) currentVizu.loadShakerMapping(EXAMPLE_MAPPING);
        });
    };

    // Remove preset left-click centering behaviour (added in NGL v2.0.0-dev.11).
    stage.mouseControls.remove("clickPick-left");

    let buttons = document.getElementsByClassName("new-bead");
    for (const button of buttons) button.disabled = true;

    const clearDialog = document.getElementById('clear-beads-dialog');
    document.getElementById('clear-beads-btn').onclick    = () => clearDialog.showModal();
    document.getElementById('clear-beads-close').onclick  = () => clearDialog.close();
    document.getElementById('clear-beads-cancel').onclick = () => clearDialog.close();
    document.getElementById('clear-beads-confirm').onclick = () => {
        if (currentVizu) {
            currentVizu.collection.clearBeads();
            currentVizu.collection.newBead();
            currentVizu.updateSelection();
        }
        clearDialog.close();
    };
    clearDialog.addEventListener('click', e => { if (e.target === clearDialog) clearDialog.close(); });

    const pasteDialog = document.getElementById('paste-mapping-dialog');
    const pasteArea   = document.getElementById('mapping-paste-area');

    document.getElementById('load-mapping-btn').onclick = () => {
        pasteArea.value = '';
        pasteDialog.showModal();
        pasteArea.focus();
    };
    document.getElementById('paste-dialog-close').onclick  = () => pasteDialog.close();
    document.getElementById('paste-dialog-cancel').onclick = () => pasteDialog.close();
    document.getElementById('paste-dialog-apply').onclick  = () => {
        const text = pasteArea.value.trim();
        if (text && currentVizu) currentVizu.loadShakerMapping(text);
        pasteDialog.close();
    };
    pasteDialog.addEventListener('click', e => { if (e.target === pasteDialog) pasteDialog.close(); });

    document.getElementById('recenter').onclick = () => {
        stage.setParameters({ clipNear: 0, clipFar: 100, fogNear: 50, fogFar: 100 });
        stage.autoView();
    };

    document.getElementById('save-image').onclick = () => {
        stage.makeImage({ factor: 2, antialias: true, trim: false, transparent: false })
            .then(blob => {
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'cgbuilder.png';
                a.click();
                URL.revokeObjectURL(url);
            });
    };

    initTheme(stage);
    initTabs();
    initNavbar(stage);
}

window.onload = main;
