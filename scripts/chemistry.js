const _ATOM_VALENCE = {C:4,N:3,O:2,S:2,P:3,F:1,CL:1,BR:1,I:1,H:1,SE:2,SI:4,B:3};
const _ORGANIC = new Set(['B','C','N','O','P','S','F','CL','BR','I']);

// Reference bond lengths (Å) as [single, double, triple]; null = not applicable.
// Used to infer bond order from 3D geometry when the input file carries no
// bond-order metadata (the common case for PDB/GRO). Keyed by element pair
// sorted alphabetically.
const _BOND_REF = {
    'C-C': [1.54, 1.34, 1.20],
    'C-N': [1.47, 1.29, 1.16],
    'C-O': [1.43, 1.23, null],
    'C-S': [1.82, 1.60, null],
    'C-P': [1.84, 1.66, null],
    'N-N': [1.45, 1.25, 1.10],
    'N-O': [1.41, 1.21, null],
    'N-S': [1.71, 1.54, null],
    'N-P': [1.65, 1.55, null],
    'O-O': [1.48, 1.21, null],
    'O-P': [1.63, 1.50, null],
    'O-S': [1.57, 1.44, null],
    'P-S': [2.12, 1.95, null],
};

function _geometryBondOrder(elA, elB, dist) {
    const ref = _BOND_REF[[elA, elB].sort().join('-')];
    if (!ref || !dist) return { order: 1, single: ref ? ref[0] : Infinity };
    let best = 1, bestErr = Infinity;
    for (let o = 1; o <= 3; o++) {
        const r = ref[o - 1];
        if (r == null) continue;
        const err = Math.abs(dist - r);
        if (err < bestErr) { bestErr = err; best = o; }
    }
    return { order: best, single: ref[0] };
}

// Perceive whole-molecule chemistry from geometry (and any NGL bond orders that
// happen to be present). Returns:
//   aromaticAtoms : Set of global atom indices in an aromatic ring
//   bondOrders    : Map "a-b" (a<b) -> order, for NON-aromatic multiple bonds
//
// Bond orders come from comparing each heavy-heavy bond length to reference
// single/double/triple lengths. Aromaticity = ring topology (prune + cycles)
// where every ring member is sp2 (carries a shortened/multiple bond) or a
// 2-connected heteroatom. Non-aromatic double/triple bonds are assigned greedily
// (most-shortened first, at most one multiple bond per atom) so resonance/amide
// carbons don't get an impossible double-double valence (C=O wins, C-N stays
// single). Aromatic bonds are left to the per-fragment kekuliser.
export function perceiveChemistry(structure) {
    const aromaticAtoms = new Set();
    const bondOrders = new Map();
    if (!structure || typeof structure.eachAtom !== 'function') {
        return { aromaticAtoms, bondOrders };
    }

    // Element + position per atom. NGL reuses one AtomProxy during iteration, so
    // read values into plain objects immediately rather than storing the proxy.
    const element = new Map();
    const pos = new Map();
    structure.eachAtom(a => {
        element.set(a.index, (a.element || 'C').toUpperCase());
        const v = typeof a.positionToVector3 === 'function'
            ? a.positionToVector3() : { x: a.x, y: a.y, z: a.z };
        pos.set(a.index, { x: v.x, y: v.y, z: v.z });
    });

    // Heavy-atom graph + per-bond geometry classification (each bond once).
    const adj = new Map();
    const bonds = [];          // { a, b, order, single, dist, aromaticNGL }
    const seen = new Set();
    const sp2 = new Set();
    const link = (u, v) => {
        if (!adj.has(u)) adj.set(u, new Set());
        adj.get(u).add(v);
    };
    structure.eachAtom(atom => {
        const idx = atom.index;
        if (element.get(idx) === 'H') return;
        if (atom.aromatic) sp2.add(idx);  // trust NGL's flag when present
        if (typeof atom.eachBond !== 'function') return;
        atom.eachBond(bond => {
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const other = i1 === idx ? i2 : (i2 === idx ? i1 : -1);
            if (other < 0 || element.get(other) === 'H') return;
            link(idx, other); link(other, idx);
            const a = Math.min(idx, other), b = Math.max(idx, other);
            const key = `${a}-${b}`;
            if (seen.has(key)) return;
            seen.add(key);
            const pa = pos.get(a), pb = pos.get(b);
            const dist = Math.hypot(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z);
            let { order, single } = _geometryBondOrder(element.get(a), element.get(b), dist);
            // Respect explicit NGL bond orders when the file does provide them.
            const ngl = bond.bondOrder;
            if (ngl === 2 || ngl === 3) order = Math.max(order, ngl);
            const aromaticNGL = ngl === 1.5 || (ngl >= 4 && ngl < 100);
            if (order >= 2 || aromaticNGL) { sp2.add(a); sp2.add(b); }
            bonds.push({ a, b, order, single, dist, aromaticNGL });
        });
    });

    // Prune terminal chains: atoms with degree <= 1 cannot be in a ring.
    const degree = new Map();
    for (const [k, s] of adj) degree.set(k, s.size);
    const removed = new Set();
    let changed = true;
    while (changed) {
        changed = false;
        for (const [k, d] of degree) {
            if (removed.has(k) || d > 1) continue;
            removed.add(k);
            changed = true;
            for (const nb of adj.get(k)) {
                if (!removed.has(nb)) degree.set(nb, degree.get(nb) - 1);
            }
        }
    }
    const ringAdj = new Map();
    for (const k of adj.keys()) {
        if (removed.has(k)) continue;
        ringAdj.set(k, [...adj.get(k)].filter(n => !removed.has(n)));
    }

    // Fundamental cycles via BFS spanning tree (handles fused rings correctly).
    // DFS back-edge detection misses cycles in polycyclic systems where the
    // closing atom of a secondary ring is already marked "done" — the BFS
    // approach finds every ring by treating each non-tree edge as a cycle.
    const cycles = [];
    {
        const bfsParent = new Map();
        const bfsQ = [];
        for (const start of ringAdj.keys()) {
            if (bfsParent.has(start)) continue;
            bfsParent.set(start, -1);
            bfsQ.push(start);
            let qi = 0;
            while (qi < bfsQ.length) {
                const u = bfsQ[qi++];
                for (const v of ringAdj.get(u)) {
                    if (!bfsParent.has(v)) { bfsParent.set(v, u); bfsQ.push(v); }
                }
            }
        }
        const treeSet = new Set();
        for (const [n, p] of bfsParent) {
            if (p >= 0) treeSet.add(`${Math.min(n, p)}-${Math.max(n, p)}`);
        }
        const seenEdge = new Set();
        for (const u of ringAdj.keys()) {
            for (const v of ringAdj.get(u)) {
                if (v <= u) continue;
                const ekey = `${u}-${v}`;
                if (treeSet.has(ekey) || seenEdge.has(ekey)) continue;
                seenEdge.add(ekey);
                // Non-tree edge u-v: find paths to LCA in spanning tree.
                const pathU = [], pathV = [];
                for (let x = u; x >= 0; x = bfsParent.get(x) ?? -1) pathU.push(x);
                for (let x = v; x >= 0; x = bfsParent.get(x) ?? -1) pathV.push(x);
                const setU = new Map(pathU.map((n, i) => [n, i]));
                let lca = -1, vIdx = -1;
                for (let i = 0; i < pathV.length; i++) {
                    if (setU.has(pathV[i])) { lca = pathV[i]; vIdx = i; break; }
                }
                if (lca < 0) continue;
                const cycle = [
                    ...pathU.slice(0, setU.get(lca) + 1),
                    ...pathV.slice(0, vIdx).reverse(),
                ];
                if (cycle.length >= 3) cycles.push(cycle);
            }
        }
    }

    const isHetero = idx => ['N', 'O', 'S'].includes(element.get(idx));
    for (const cyc of cycles) {
        if (cyc.length < 4 || cyc.length > 7) continue;
        if (cyc.every(idx => sp2.has(idx) || isHetero(idx))) {
            for (const idx of cyc) aromaticAtoms.add(idx);
        }
    }

    // Assign non-aromatic multiple bonds greedily: most-shortened first, at most
    // one per atom. Aromatic bonds are excluded — kekulisation handles those.
    // The one-per-atom cap stops resonance/amide carbons getting an impossible
    // double-double valence: e.g. for an amide the C=O (more shortened) wins and
    // the C–N is left single.
    const ranked = bonds
        .filter(bd => bd.order >= 2 && !(aromaticAtoms.has(bd.a) && aromaticAtoms.has(bd.b)))
        .sort((x, y) => (y.single - y.dist) - (x.single - x.dist));
    const usedMul = new Set();
    for (const bd of ranked) {
        if (usedMul.has(bd.a) || usedMul.has(bd.b)) continue;
        bondOrders.set(`${bd.a}-${bd.b}`, bd.order);
        usedMul.add(bd.a); usedMul.add(bd.b);
    }
    return { aromaticAtoms, bondOrders };
}

// aromaticNotation=true: write aromatic atoms as lowercase (c, n, o…), skip
// kekulization. Used to generate AutoMartini-style lookup keys ("cc", "ccc"…)
// that match open-chain aromatic entries in the table. Never feed this output
// to RDKit — it rejects open-chain aromatic SMILES.
//
// aromaticAtoms: optional Set of global atom indices known to be aromatic (from
// perceiveChemistry). When provided it overrides the local flag/bond-order
// heuristic, which is unreliable for Kekulé-form input.
// bondOrders: optional Map "a-b"(a<b) -> order for non-aromatic multiple bonds
// (from perceiveChemistry's geometry inference). Overrides NGL bond orders so
// carbonyls/imines/alkenes render correctly even from PDB/GRO without bond data.
export function fragmentToSmiles(beadAtoms,
        { aromaticNotation = false, aromaticAtoms = null, bondOrders = null,
          startIndex = null } = {}) {
    const heavy = beadAtoms.filter(a => (a.element || 'C').toUpperCase() !== 'H');
    if (heavy.length === 0) return null;

    const atomSet = new Set(heavy.map(a => a.index));

    // Determine which bead atoms are aromatic. Prefer the caller-supplied set
    // (whole-molecule ring perception); otherwise fall back to the local
    // flag/bond-order heuristic (works only when NGL encodes aromaticity).
    const aromaticSet = new Set();
    if (aromaticAtoms) {
        for (const atom of heavy) {
            if (aromaticAtoms.has(atom.index)) aromaticSet.add(atom.index);
        }
    } else {
        for (const atom of heavy) {
            if (atom.aromatic) { aromaticSet.add(atom.index); continue; }
            if (typeof atom.eachBond !== 'function') continue;
            let found = false;
            atom.eachBond(b => {
                const bo = b.bondOrder;
                if (bo === 1.5 || (bo >= 4 && bo < 100)) found = true;
            });
            if (found) aromaticSet.add(atom.index);
        }
    }

    // Collect each internal bond once (undirected), flagging aromatic bonds.
    // NGL may encode aromaticity as bond order 4 or 1.5, or only via the per-atom
    // aromatic flag, so check both.
    const edges = new Map(); // "min-max" -> { a, b, order, aromatic }
    for (const atom of heavy) {
        if (typeof atom.eachBond !== 'function') continue;
        atom.eachBond(bond => {
            // Identify the other endpoint robustly: NGL may order atomIndex1/2
            // with the lower index first (not necessarily the current atom first).
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const otherIdx = i1 === atom.index ? i2 : (i2 === atom.index ? i1 : -1);
            if (otherIdx < 0 || !atomSet.has(otherIdx)) return;
            const a = Math.min(atom.index, otherIdx);
            const b = Math.max(atom.index, otherIdx);
            const key = `${a}-${b}`;
            if (edges.has(key)) return;
            const raw = bond.bondOrder || 1;
            const aromatic = raw >= 4 || raw === 1.5
                || (aromaticSet.has(a) && aromaticSet.has(b));
            // Bond order: aromatic bonds get kekulised below; otherwise prefer
            // the geometry-inferred order, falling back to NGL's bond order.
            let order;
            if (aromatic) order = 1;
            else if (bondOrders && bondOrders.has(key)) order = bondOrders.get(key);
            else order = Math.round(raw);
            edges.set(key, { a, b, order, aromatic });
        });
    }

    // Kekulise aromatic bonds (Kekulé mode only): a greedy maximum matching
    // turns alternating aromatic bonds into double bonds so sp2 atoms keep
    // their unsaturation. Skipped in aromaticNotation mode where aromatic atoms
    // are written lowercase and bonds are implicit (order stays 1).
    if (!aromaticNotation) {
        const doubled = new Set();
        for (const edge of edges.values()) {
            if (edge.aromatic && !doubled.has(edge.a) && !doubled.has(edge.b)) {
                edge.order = 2;
                doubled.add(edge.a);
                doubled.add(edge.b);
            }
        }
    }

    // Per-atom data: internal bonds and implicit H count. Any valence not used by
    // an internal bond is capped with hydrogen (external heavy bonds and original
    // explicit hydrogens alike).
    const data = new Map();
    for (const atom of heavy) {
        data.set(atom.index, {
            el: (atom.element || 'C').toUpperCase(),
            charge: atom.formalCharge ?? 0,
            hCount: 0,
            internalBonds: [],
        });
    }
    const intSum = new Map();
    for (const edge of edges.values()) {
        data.get(edge.a).internalBonds.push({ toIdx: edge.b, order: edge.order });
        data.get(edge.b).internalBonds.push({ toIdx: edge.a, order: edge.order });
        intSum.set(edge.a, (intSum.get(edge.a) || 0) + edge.order);
        intSum.set(edge.b, (intSum.get(edge.b) || 0) + edge.order);
    }
    for (const atom of heavy) {
        const d = data.get(atom.index);
        const baseVal = _ATOM_VALENCE[d.el] ?? 4;
        d.hCount = Math.max(0, baseVal + d.charge - (intSum.get(atom.index) || 0));
    }

    // DFS to find ring-closure back-edges and assign ring-closure digits
    const visited = new Set();
    const inStack = new Set();
    const seenEdges = new Set();
    const ringClosures = new Map(); // atomIdx -> [{digit, writeBond}]
    let digitCounter = 1;

    function findBackEdges(idx, parentIdx) {
        visited.add(idx); inStack.add(idx);
        for (const { toIdx, order } of (data.get(idx)?.internalBonds ?? [])) {
            const edgeKey = `${Math.min(idx, toIdx)}-${Math.max(idx, toIdx)}`;
            if (toIdx === parentIdx || seenEdges.has(edgeKey)) continue;
            seenEdges.add(edgeKey);
            if (inStack.has(toIdx)) {
                const digit = digitCounter++;
                if (!ringClosures.has(toIdx)) ringClosures.set(toIdx, []);
                if (!ringClosures.has(idx))   ringClosures.set(idx,   []);
                ringClosures.get(toIdx).push({ digit, writeBond: false });
                ringClosures.get(idx).push({ digit, order, writeBond: true });
            } else if (!visited.has(toIdx)) {
                findBackEdges(toIdx, idx);
            }
        }
        inStack.delete(idx);
    }

    const startIdx = (startIndex != null && atomSet.has(startIndex))
        ? startIndex : heavy[0].index;
    findBackEdges(startIdx, -1);
    visited.clear();

    function bondChar(order) {
        return order === 2 ? '=' : order === 3 ? '#' : '';
    }

    function atomToken(idx) {
        const d = data.get(idx);
        if (aromaticNotation && aromaticSet.has(idx)) {
            return d.el.charAt(0).toLowerCase() + d.el.slice(1).toLowerCase();
        }
        const sym = d.el === 'CL' ? 'Cl' : d.el === 'BR' ? 'Br'
                  : d.el.charAt(0) + d.el.slice(1).toLowerCase();

        // hCount is computed from the same valence table SMILES uses for implicit H,
        // so brackets are never needed just because hCount > 0 for organic atoms.
        const needsBracket = d.charge !== 0 || !_ORGANIC.has(d.el);
        if (!needsBracket) return sym;

        let inner = sym;
        if (d.hCount === 1) inner += 'H';
        else if (d.hCount > 1) inner += `H${d.hCount}`;
        if (d.charge > 0) inner += d.charge === 1 ? '+' : `+${d.charge}`;
        else if (d.charge < 0) inner += d.charge === -1 ? '-' : `${d.charge}`;
        return `[${inner}]`;
    }

    function closureSuffix(idx) {
        return (ringClosures.get(idx) ?? []).map(({ digit, order, writeBond }) => {
            const b = writeBond ? bondChar(order) : '';
            return b + (digit >= 10 ? `%${digit}` : `${digit}`);
        }).join('');
    }

    function dfs(idx) {
        visited.add(idx);
        let smi = atomToken(idx) + closureSuffix(idx);

        const children = (data.get(idx)?.internalBonds ?? [])
            .filter(({ toIdx }) => !visited.has(toIdx));

        if (children.length === 0) return smi;
        // Branches must immediately follow the parent atom. Write every child
        // except the last as a parenthesised branch, then the last child inline
        // as the main-chain continuation. Children visited during a sibling's
        // recursion (rings) are skipped; their ring-closure digit is emitted via
        // closureSuffix instead.
        for (let i = 0; i < children.length - 1; i++) {
            const { toIdx, order } = children[i];
            if (visited.has(toIdx)) continue;
            smi += `(${bondChar(order)}${dfs(toIdx)})`;
        }
        const last = children[children.length - 1];
        if (!visited.has(last.toIdx)) {
            smi += bondChar(last.order) + dfs(last.toIdx);
        }
        return smi;
    }

    return dfs(startIdx);
}

const _HETERO_VALENCE = { N: 3, O: 2, S: 2 };

// Number of hydrogens actually attached to a heteroatom in the real structure.
// hasExplicitH: whether the loaded molecule carries explicit hydrogen atoms.
export function realHydrogenCount(atom, hasExplicitH) {
    const el = (atom.element || '').toUpperCase();
    const baseVal = _HETERO_VALENCE[el];
    if (baseVal === undefined) return 0;

    const charge = atom.formalCharge ?? 0;
    const structure = atom.structure;
    let explicitH = 0;
    let heavySum = 0;
    let aromaticRingBonds = 0;

    if (typeof atom.eachBond === 'function') {
        atom.eachBond((bond) => {
            const otherIdx = bond.atomIndex1 === atom.index
                ? bond.atomIndex2 : bond.atomIndex1;
            const other = structure.getAtomProxy(otherIdx);
            if ((other.element || '').toUpperCase() === 'H') { explicitH += 1; return; }
            let order = bond.bondOrder || 1;
            // Treat as an aromatic ring bond if both atoms are flagged aromatic,
            // OR the bond order encodes aromaticity (NGL may use 4 or 1.5).
            const aromatic = (!!atom.aromatic && !!other.aromatic)
                || order >= 4 || order === 1.5;
            if (aromatic) { aromaticRingBonds += 1; order = 1; }
            heavySum += order;
        });
    }

    // If the molecule carries explicit hydrogens, the count of H neighbours is
    // definitive — even when it is zero (e.g. a pyridine-type aromatic N).
    if (hasExplicitH) return explicitH;

    // No explicit H: estimate implicit H from valence. An aromatic ring atom
    // consumes one extra valence unit beyond its (normalised) single ring bonds,
    // so pyridine N -> 0 H and benzene C -> 1 H. (A pyrrole-type N-H cannot be
    // told apart from a pyridine N without explicit hydrogens.)
    let used = heavySum;
    if (aromaticRingBonds >= 2) used += 1;
    return Math.max(0, baseVal + charge - used);
}

export function beadDonorCount(bead, hasExplicitH) {
    let donors = 0;
    for (const atom of bead.atoms) {
        const el = (atom.element || '').toUpperCase();
        if ((el === 'N' || el === 'O' || el === 'S')
            && realHydrogenCount(atom, hasExplicitH) > 0) {
            donors += 1;
        }
    }
    return donors;
}

export function structureHasHydrogens(structure) {
    let found = false;
    if (structure && typeof structure.eachAtom === 'function') {
        structure.eachAtom((atom) => {
            if ((atom.element || '').toUpperCase() === 'H') found = true;
        });
    }
    return found;
}
