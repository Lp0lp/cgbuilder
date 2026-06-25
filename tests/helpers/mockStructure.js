// Minimal duck-typed stand-ins for NGL's Structure/AtomProxy/Bond interfaces,
// just enough surface for chemistry.js/sasa.js/fileformats.js to operate on.
export function buildStructure(atomDefs, bondDefs = []) {
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

    const bondsByAtom = new Map(atoms.map((a) => [a.index, []]));
    for (const bd of bondDefs) {
        const bond = { atomIndex1: bd.a, atomIndex2: bd.b, bondOrder: bd.order ?? 1 };
        bondsByAtom.get(bd.a).push(bond);
        bondsByAtom.get(bd.b).push(bond);
    }

    const structure = {};
    const proxies = atoms.map((atom) => ({
        ...atom,
        structure,
        positionToVector3: () => ({ x: atom.x, y: atom.y, z: atom.z }),
        eachBond: (cb) => { for (const b of bondsByAtom.get(atom.index)) cb(b); },
    }));

    structure.eachAtom = (cb) => { for (const p of proxies) cb(p); };
    structure.getAtomProxy = (idx) => proxies[idx];
    return structure;
}

// Atoms of a regular N-gon (side length = bondLen) centred at (cx, cy), in the
// XY plane, returned as atomDefs ready for buildStructure. Used to build
// aromatic test rings without hand-computing trigonometry per test.
export function regularPolygon(n, bondLen, element, { cx = 0, cy = 0, startAngleDeg = 0 } = {}) {
    const R = bondLen / (2 * Math.sin(Math.PI / n));
    const pts = [];
    for (let k = 0; k < n; k++) {
        const theta = (startAngleDeg * Math.PI / 180) + (2 * Math.PI * k) / n;
        pts.push({ element, x: cx + R * Math.cos(theta), y: cy + R * Math.sin(theta), z: 0 });
    }
    return pts;
}

export function ringBonds(indices, order) {
    const bonds = [];
    for (let i = 0; i < indices.length; i++) {
        bonds.push({ a: indices[i], b: indices[(i + 1) % indices.length], order });
    }
    return bonds;
}
