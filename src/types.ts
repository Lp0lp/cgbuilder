/* ===========================================================================
   Shared type definitions
   ===========================================================================
   Minimal, duck-typed interfaces for the third-party objects this app talks
   to (NGL's structure/stage/shape API, RDKit_minimal's molecule wrapper) and
   for the plain data shapes passed between our own modules (perceived
   chemistry, parsed mapping definitions, bead-type prediction inputs). Only
   the members the code actually uses are declared — these are deliberately
   partial views of much larger real APIs, not full type stubs. */

/** An NGL.Vector3 — the subset of its API used for weighted-centre math. */
export interface Vec3 {
    x: number;
    y: number;
    z: number;
    add(v: Vec3): this;
    divideScalar(s: number): this;
}

/** One bond, as exposed by NGL's BondProxy (order is directly settable). */
export interface BondProxy {
    atomIndex1: number;
    atomIndex2: number;
    bondOrder: number;
}

/** One atom, as exposed by NGL's AtomProxy (plus our mock in tests). */
export interface AtomProxy {
    index: number;
    element?: string;
    atomname: string;
    resname?: string;
    resno: number;
    chainname?: string;
    x: number;
    y: number;
    z: number;
    structure: Structure;
    positionToVector3(): Vec3;
    eachBond(cb: (bond: BondProxy) => void): void;
}

/** An NGL Structure — the iteration/lookup surface our modules rely on. */
export interface Structure {
    eachAtom(cb: (atom: AtomProxy) => void): void;
    getAtomProxy(index: number): AtomProxy;
}

/** NGL representation parameters — an untyped bag forwarded straight to NGL. */
export type RepresentationParams = Record<string, unknown>;

/** A single NGL representation handle (ball+stick, surface, label, buffer). */
export interface NglRepresentation {
    setSelection(selection: string): void;
    setVisibility(visible: boolean): void;
}

/** Base NGL component: anything that can carry representations. */
export interface Component {
    addRepresentation(type: string, params?: RepresentationParams): NglRepresentation;
}

/** An NGL StructureComponent: a component backed by a real Structure. */
export interface StructureComponent extends Component {
    structure: Structure;
    autoView(): void;
    rebuildRepresentations(): void;
}

/** NGL PickingProxy for a viewer click (atom is absent for empty-space clicks). */
export interface PickingProxy {
    atom?: AtomProxy;
    mouse?: { shiftKey: boolean };
}

/** An NGL.Shape under construction (CG bead spheres + labels). */
export interface Shape {
    addSphere(position: Vec3 | number[], color: number[], radius: number, name?: string): void;
    addText(position: number[], color: number[], size: number, text: string): void;
}

/** The NGL Stage — viewer lifecycle, component management, IO, image export. */
export interface Stage {
    signals: {
        clicked: {
            removeAll(): void;
            add(cb: (proxy: PickingProxy | undefined) => void): void;
        };
    };
    mouseControls: { remove(name: string): void };
    removeAllComponents(): void;
    removeComponent(component: Component): void;
    addComponentFromObject(shape: Shape): Component;
    loadFile(file: File | Blob, params?: object): Promise<StructureComponent>;
    setParameters(params: object): void;
    autoView(): void;
    handleResize(): void;
    makeImage(params: object): Promise<Blob>;
}

/** The `NGL` browser global (loaded via `<script>` tag in index.html). */
export interface NGLStatic {
    Vector3: new (x?: number, y?: number, z?: number) => Vec3;
    Stage: new (idOrElement: string | HTMLElement) => Stage;
    Shape: new (name: string, params?: object) => Shape;
}

/** The initialized RDKit_minimal module (re-exported from the npm package
 *  so consumers keep a single `./types.js` import surface; see rdkit.ts). */
export type { RDKitModule } from '@rdkit/rdkit';

/**
 * Whole-molecule chemistry perceived from connectivity + explicit hydrogens
 * (see chemistry.ts's perceiveChemistry). Consumers must check `available`
 * before trusting any other field.
 */
export interface Chemistry {
    available: boolean;
    ringAtoms: Set<number>;
    aromaticAtoms: Set<number>;
    branchAtoms: Set<number>;
    bondOrders: Map<string, number>;
    charges: Map<number, number>;
    hNeighbors: Map<number, number>;
}

/** One bead definition parsed from a Shaker mapping (see parseShakerMapping). */
export interface BeadDef {
    name: string;
    type: string;
    charge: number;
    atoms: string[];
}

/** Inputs to determineBeadType (see prediction.ts). */
export interface BeadTypeProps {
    deltaF: number;
    charge: number;
    hDonors: number;
    hAcceptors: number;
    hasHalogen: boolean;
    inRing: boolean;
    weightedHeavyCount: number;
    ringOrBranched: boolean;
}
