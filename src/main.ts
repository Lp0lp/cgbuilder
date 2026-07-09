import { BeadCollection } from './bead.js';
import { Visualization } from './visualization.js';
import { readOriginalAtomNames, bondAwareRepresentationParams } from './fileformats.js';
import { EXAMPLE_PDB, EXAMPLE_MAPPING } from './example.js';
import { byId } from './dom.js';
import { NGL } from './ngl.js';
import { marked } from 'marked';
import type { Stage } from './types.js';

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
let currentVizu: Visualization | null = null;

/**
 * Load a molecule from a File into the viewer: tears down any previously
 * loaded component, creates a fresh BeadCollection + Visualization for it,
 * and loads both the structure (via NGL) and its original atom names (see
 * fileformats.ts's readOriginalAtomNames) in parallel before wiring up the
 * rest of the per-molecule UI (new-bead buttons, the 3D click handler).
 * Every one of those controls stays disabled/unregistered until the load
 * actually finishes — they all end up calling Visualization.updateSelection,
 * which dereferences `representation`, only set once attachRepresentation
 * runs below; wiring them any earlier would let a click during the load
 * (e.g. on an already-disabled-looking button, or on empty viewport space
 * before the structure has rendered) throw on a null representation.
 * @param file
 * @param stage - NGL Stage
 * @returns resolves once the molecule is fully loaded and the UI is ready —
 *   callers that need to act afterward (e.g. applying a mapping right after
 *   loading the bundled example) chain onto this
 */
function loadMoleculeFromFile(file: File, stage: Stage): Promise<void> {
    stage.removeAllComponents();
    stage.signals.clicked.removeAll();

    const collection = new BeadCollection();
    const vizu = new Visualization(collection, stage);
    currentVizu = vizu;

    const newBeadButtons = document.getElementsByClassName("new-bead") as HTMLCollectionOf<HTMLButtonElement>;
    for (const button of newBeadButtons) button.disabled = true;
    byId<HTMLButtonElement>('load-mapping-btn').disabled = true;
    byId<HTMLButtonElement>('clear-beads-btn').disabled = true;

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
            byId<HTMLButtonElement>('load-mapping-btn').disabled = false;
            byId<HTMLButtonElement>('clear-beads-btn').disabled = false;
            stage.signals.clicked.add((pickingProxy) => vizu.onClick(pickingProxy));
        })
        .catch((err) => {
            console.error("Error loading molecule or reading original atom names:", err);
        });

    return ready;
}

/**
 * The "Choose File" <input type="file">'s change handler.
 * @param event
 * @param stage - NGL Stage
 */
function loadMolecule(event: Event, stage: Stage): void {
    const files = (event.target as HTMLInputElement).files;
    if (files && files[0]) loadMoleculeFromFile(files[0], stage);
}

/** Wire up the output-tab bar (Shaker/.gro/.ndx/.map/AA SMILES) switching. */
function initTabs(): void {
    const btns = document.querySelectorAll<HTMLElement>('.tab-btn');
    const panels = document.querySelectorAll<HTMLElement>('.tab-panel');
    btns.forEach(btn => {
        btn.onclick = () => {
            btns.forEach(b => b.classList.remove('active'));
            panels.forEach(p => p.classList.remove('active'));
            btn.classList.add('active');
            byId('tab-' + btn.dataset.tab).classList.add('active');
        };
    });
}

/**
 * Wire up the dark/light theme toggle and keep the NGL stage's background
 * colour in sync with it (and vice versa, via the separate "Light BG"
 * checkbox — the two are independent: page theme is persisted to
 * localStorage and restored on load via the inline <head> script; the
 * viewer background is a per-session NGL parameter, not persisted).
 * @param stage - NGL Stage
 */
function initTheme(stage: Stage): void {
    const themeBtn = byId('theme-toggle');
    const bgToggle = byId<HTMLInputElement>('toggle-bg');

    // Apply isDark to both the NGL viewer background and the "Light BG"
    // checkbox's own state, so the two stay consistent with each other.
    function syncBg(isDark: boolean): void {
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
 * @param pageEl - the page's container, holding a
 *   `.doc-content[data-md]` element naming its markdown source file
 */
function loadDocPage(pageEl: Element): void {
    const container = pageEl.querySelector<HTMLElement>('.doc-content[data-md]');
    if (!container || container.dataset.loaded) return;
    fetch(container.dataset.md ?? '')
        .then(r => r.text())
        .then(md => { container.innerHTML = marked.parse(md, { async: false }); container.dataset.loaded = '1'; })
        .catch(() => { container.innerHTML = '<p>Could not load documentation.</p>'; });
}

/**
 * Wire up the top-level navbar (App / How-to / Guidelines / ...): switches
 * the visible `.page`, and either resizes the NGL stage (returning to the
 * App page, whose canvas may be stale after being hidden via CSS while
 * another page was active) or lazily loads that page's doc content
 * (loadDocPage).
 * @param stage - NGL Stage
 */
function initNavbar(stage: Stage): void {
    const tabs = document.querySelectorAll<HTMLElement>('.navbar-tab');
    const pages = document.querySelectorAll<HTMLElement>('.page');
    tabs.forEach(tab => {
        tab.onclick = () => {
            tabs.forEach(t => t.classList.remove('active'));
            pages.forEach(p => p.classList.remove('active'));
            tab.classList.add('active');
            const pageEl = byId('page-' + tab.dataset.page);
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
function main(): void {
    // Capture wheel events within the viewer so the page doesn't scroll when zooming.
    // https://github.com/nglviewer/ngl/issues/878#issuecomment-913504711
    const stageContainer = byId('viewport');
    window.addEventListener('wheel', (event) => {
        if (stageContainer.contains(event.target as Node)) event.preventDefault();
    }, { passive: false });

    const stage = new NGL.Stage("viewport");

    window.addEventListener("resize", () => stage.handleResize(), false);

    const mol_select = byId<HTMLInputElement>("mol-select");
    mol_select.onchange = (event) => loadMolecule(event, stage);

    byId<HTMLButtonElement>('load-example-btn').onclick = () => {
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

    const buttons = document.getElementsByClassName("new-bead") as HTMLCollectionOf<HTMLButtonElement>;
    for (const button of buttons) button.disabled = true;

    const clearDialog = byId<HTMLDialogElement>('clear-beads-dialog');
    byId<HTMLButtonElement>('clear-beads-btn').onclick    = () => clearDialog.showModal();
    byId<HTMLButtonElement>('clear-beads-close').onclick  = () => clearDialog.close();
    byId<HTMLButtonElement>('clear-beads-cancel').onclick = () => clearDialog.close();
    byId<HTMLButtonElement>('clear-beads-confirm').onclick = () => {
        if (currentVizu) {
            currentVizu.collection.clearBeads();
            currentVizu.collection.newBead();
            currentVizu.updateSelection();
        }
        clearDialog.close();
    };
    clearDialog.addEventListener('click', e => { if (e.target === clearDialog) clearDialog.close(); });

    const pasteDialog = byId<HTMLDialogElement>('paste-mapping-dialog');
    const pasteArea   = byId<HTMLTextAreaElement>('mapping-paste-area');

    byId<HTMLButtonElement>('load-mapping-btn').onclick = () => {
        pasteArea.value = '';
        pasteDialog.showModal();
        pasteArea.focus();
    };
    byId<HTMLButtonElement>('paste-dialog-close').onclick  = () => pasteDialog.close();
    byId<HTMLButtonElement>('paste-dialog-cancel').onclick = () => pasteDialog.close();
    byId<HTMLButtonElement>('paste-dialog-apply').onclick  = () => {
        const text = pasteArea.value.trim();
        if (text && currentVizu) currentVizu.loadShakerMapping(text);
        pasteDialog.close();
    };
    pasteDialog.addEventListener('click', e => { if (e.target === pasteDialog) pasteDialog.close(); });

    byId<HTMLButtonElement>('recenter').onclick = () => {
        stage.setParameters({ clipNear: 0, clipFar: 100, fogNear: 50, fogFar: 100 });
        stage.autoView();
    };

    byId<HTMLButtonElement>('save-image').onclick = () => {
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
