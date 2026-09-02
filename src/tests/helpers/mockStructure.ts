// Minimal duck-typed stand-ins for NGL's Structure/AtomProxy/Bond interfaces,
// just enough surface for chemistry.ts/sasa.ts/fileformats.ts to operate on.
import type { AtomProxy, BondProxy, Structure, Vec3 } from '../../types.js';

/** Input shape for one atom (all but `element` optional, defaulted below). */
export interface AtomDef {
    element: string;
    x?: number;
    y?: number;
    z?: number;
    formalCharge?: number;
    aromatic?: boolean;
    resname?: string;
    resno?: number;
    atomname?: string;
}

/** Input shape for one bond (order defaults to single). */
export interface BondDef {
    a: number;
    b: number;
    order?: number;
}

export function buildStructure(atomDefs: AtomDef[], bondDefs: BondDef[] = []): Structure {
    const atoms = atomDefs.map((def, index) => ({
        index,
        element: def.element,
        x: def.x ?? 0,
        y: def.y ?? 0,
        z: def.z ?? 0,
        formalCharge: def.formalCharge ?? 0,
        aromatic: def.aromatic ?? false,
        resname: def.resname ?? 'UNK',
        resno: def.resno ?? 1,
        atomname: def.atomname ?? `A${index}`,
    }));

    const bondsByAtom = new Map<number, BondProxy[]>(
        atoms.map((a): [number, BondProxy[]] => [a.index, []]),
    );
    for (const bd of bondDefs) {
        const bond: BondProxy = { atomIndex1: bd.a, atomIndex2: bd.b, bondOrder: bd.order ?? 1 };
        bondsByAtom.get(bd.a)!.push(bond);
        bondsByAtom.get(bd.b)!.push(bond);
    }

    const structure = {} as Structure;
    const proxies: AtomProxy[] = atoms.map((atom) => ({
        ...atom,
        structure,
        positionToVector3: (): Vec3 => ({ x: atom.x, y: atom.y, z: atom.z }) as Vec3,
        eachBond: (cb: (bond: BondProxy) => void) => { for (const b of bondsByAtom.get(atom.index)!) cb(b); },
    }));

    structure.eachAtom = (cb: (atom: AtomProxy) => void) => { for (const p of proxies) cb(p); };
    structure.getAtomProxy = (idx: number) => proxies[idx];
    return structure;
}

// Atoms of a regular N-gon (side length = bondLen) centred at (cx, cy), in the
// XY plane, returned as atomDefs ready for buildStructure. Used to build
// aromatic test rings without hand-computing trigonometry per test.
export function regularPolygon(
    n: number,
    bondLen: number,
    element: string,
    { cx = 0, cy = 0, startAngleDeg = 0 }: { cx?: number; cy?: number; startAngleDeg?: number } = {},
): AtomDef[] {
    const R = bondLen / (2 * Math.sin(Math.PI / n));
    const pts: AtomDef[] = [];
    for (let k = 0; k < n; k++) {
        const theta = (startAngleDeg * Math.PI / 180) + (2 * Math.PI * k) / n;
        pts.push({ element, x: cx + R * Math.cos(theta), y: cy + R * Math.sin(theta), z: 0 });
    }
    return pts;
}

export function ringBonds(indices: number[], order?: number): BondDef[] {
    const bonds: BondDef[] = [];
    for (let i = 0; i < indices.length; i++) {
        bonds.push({ a: indices[i], b: indices[(i + 1) % indices.length], order });
    }
    return bonds;
}
