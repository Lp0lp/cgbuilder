import { BeadCollection } from './bead.js';
import { Visualization } from './visualization.js';
import { readOriginalAtomNames, bondAwareRepresentationParams } from './fileformats.js';
import { EXAMPLE_PDB, EXAMPLE_MAPPING } from './example.js';

/* ===========================================================================
   main — entry point and page-level wiring
   ===========================================================================
   Bootstraps the NGL stage, creates a fresh BeadCollection + Visualization
   per loaded molecule, and wires up everything that isn't specific to one
   molecule (theme, the output-tab bar, the top-level navbar/doc pages, the
   clear-beads and paste-mapping dialogs). Pure UI wiring — no algorithm
   lives here... */

// The active Visualization instance, since the app only ever has one
// molecule loaded at a time. Module-level rather than passed around, so
// page-level handlers set up once in main() (the clear-beads/paste-mapping
// dialogs) can reach whichever molecule is currently loaded without each
// needing their own reference threaded through.
let currentVizu = null;

/**
 * Load a molecule from a File into the viewer: tears down any previously
 * loaded component, creates a fresh BeadCollection + Visualization for it,
 * and loads both the structure (via NGL) and its original atom names (see
 * fileformats.js's readOriginalAtomNames) in parallel before wiring up the
 * rest of the per-molecule UI (new-bead buttons, the 3D click handler).
 * Every one of those controls stays disabled/unregistered until the load
 * actually finishes — they all end up calling Visualization.updateSelection,
 * which dereferences `representation`, only set once attachRepresentation
 * runs below; wiring them any earlier would let a click during the load
 * (e.g. on an already-disabled-looking button, or on empty viewport space
 * before the structure has rendered) throw on a null representation.
 * @param {File} file
 * @param {object} stage - NGL Stage
 * @returns {Promise<void>} resolves once the molecule is fully loaded and
 *   the UI is ready — callers that need to act afterward (e.g. applying a
 *   mapping right after loading the bundled example) chain onto this
 */
function loadMoleculeFromFile(file, stage) {
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();

    let collection = new BeadCollection();
    let vizu = new Visualization(collection, stage);
    currentVizu = vizu;

    let newBeadButtons = document.getElementsByClassName("new-bead");
    for (const button of newBeadButtons) button.disabled = true;
    document.getElementById('load-mapping-btn').disabled = true;
    document.getElementById('clear-beads-btn').disabled = true;

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

            for (const button of newBeadButtons) {
                button.onclick = (event) => vizu.onNewBead(event);
                button.disabled = false;
            }
            document.getElementById('load-mapping-btn').disabled = false;
            document.getElementById('clear-beads-btn').disabled = false;
            stage.signals.clicked.add((pickingProxy) => vizu.onClick(pickingProxy));
        })
        .catch((err) => {
            console.error("Error loading molecule or reading original atom names:", err);
        });

    return ready;
}

/**
 * The "Choose File" <input type="file">'s change handler.
 * @param {Event} event
 * @param {object} stage - NGL Stage
 */
function loadMolecule(event, stage) {
    loadMoleculeFromFile(event.target.files[0], stage);
}

/** Wire up the output-tab bar (Shaker/.gro/.ndx/.map/AA SMILES) switching. */
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

/**
 * Wire up the dark/light theme toggle and keep the NGL stage's background
 * colour in sync with it (and vice versa, via the separate "Light BG"
 * checkbox — the two are independent: page theme is persisted to
 * localStorage and restored on load via the inline <head> script; the
 * viewer background is a per-session NGL parameter, not persisted).
 * @param {object} stage - NGL Stage
 */
function initTheme(stage) {
    const themeBtn = document.getElementById('theme-toggle');
    const bgToggle = document.getElementById('toggle-bg');

    // Apply isDark to both the NGL viewer background and the "Light BG"
    // checkbox's own state, so the two stay consistent with each other.
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

/**
 * Lazily fetch and render one navbar page's markdown content (How-to,
 * Guidelines, ...) into HTML via `marked`, the first time that page is
 * visited. `container.dataset.loaded` memoizes this so revisiting a page
 * doesn't re-fetch/re-render it.
 * @param {Element} pageEl - the page's container, holding a
 *   `.doc-content[data-md]` element naming its markdown source file
 */
function loadDocPage(pageEl) {
    const container = pageEl.querySelector('.doc-content[data-md]');
    if (!container || container.dataset.loaded) return;
    fetch(container.dataset.md)
        .then(r => r.text())
        .then(md => { container.innerHTML = marked.parse(md); container.dataset.loaded = '1'; })
        .catch(() => { container.innerHTML = '<p>Could not load documentation.</p>'; });
}

/**
 * Wire up the top-level navbar (App / How-to / Guidelines / ...): switches
 * the visible `.page`, and either resizes the NGL stage (returning to the
 * App page, whose canvas may be stale after being hidden via CSS while
 * another page was active) or lazily loads that page's doc content
 * (loadDocPage).
 * @param {object} stage - NGL Stage
 */
function initNavbar(stage) {
    const tabs = document.querySelectorAll('.navbar-tab');
    const pages = document.querySelectorAll('.page');
    tabs.forEach(tab => {
        tab.onclick = () => {
            tabs.forEach(t => t.classList.remove('active'));
            pages.forEach(p => p.classList.remove('active'));
            tab.classList.add('active');
            const pageEl = document.getElementById('page-' + tab.dataset.page);
            pageEl.classList.add('active');
            if (tab.dataset.page === 'app') stage.handleResize();
            else loadDocPage(pageEl);
        };
    });
}

/**
 * App entry point (see window.onload below): creates the NGL stage, wires
 * up every control that doesn't depend on a molecule being loaded yet (file
 * input, the bundled-example button, the clear-beads and paste-mapping
 * dialogs, recenter/save-image), and delegates the rest to initTheme/
 * initTabs/initNavbar.
 */
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
    // Remove NGL's built-in ctrl-click distance/angle/dihedral measurement
    // tool — it fires alongside our own atom-click-to-bead handling and
    // looks like unrelated, confusing behaviour in this app.
    stage.mouseControls.remove("clickPick-ctrl-left");

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

// Run once the whole page (not just the DOM) has loaded, since main() reads
// layout-dependent state (NGL.Stage sizes itself from the viewport element).
window.onload = main;
