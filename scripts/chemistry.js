/* ===========================================================================
   Chemistry perception — bond order, aromaticity, formal charge
   ===========================================================================
   ALGORITHM. PDB/GRO files carry atomic connectivity but no bond *order* —
   the file never says whether a bond is single, double, or triple. This
   module used to guess bond order from 3D bond length (a shortened bond
   "looks like" a double bond). That guess is unreliable for ring systems:
   resonance blurs aromatic bond lengths, and it had no concept of "this
   heteroatom's lone pair completes the ring, it can never take a double
   bond" — which produced real, confirmed bugs (a furan-type oxygen could
   get assigned a double bond, making the fragment invalid).

   This rewrite instead tries to follow the approach from xyz2mol (Kim & Kim,
   2015, https://doi.org/10.1002/bkcs.10334; reference implementation
   https://github.com/jensengroup/xyz2mol): derive bond order and formal
   charge from valence bookkeeping alone, given a structure with EXPLICIT
   hydrogen atoms (see "why explicit hydrogens are required" below).

     1. Every atom has a target valence (C=4, N=3, O=2, S=2, P=3, ... — see
        _VALENCE_OPTIONS, which also lists higher options for hypervalent
        atoms like phosphate P or sulfonate S).
     2. "Deficit" = target valence − current bond-order sum (every bond,
        including to real H atoms, starts at order 1). An atom with
        deficit > 0 has room for one more bond.
     3. A maximum matching (graph theory: the largest set of bonds with no
        shared endpoint) is found among every bonded pair of deficit>0
        atoms. Each matched bond becomes a double bond. This picks ALL the
        double bonds in a ring simultaneously and consistently, rather than
        guessing bond-by-bond from geometry. A furan/thiophene ring
        oxygen/sulfur (two single ring bonds, no H) already has zero
        deficit and is automatically excluded from the matching — not
        through a hard-coded "exclude O/S" rule, simply because it has
        nothing left to fill.
     4. Triple bonds: a matched pair that still has exactly 1 unit of
        deficit left on both sides, with no other eligible neighbour for
        either atom, escalates to a triple bond (nitriles, alkynes).
     5. Formal charge falls out of the same bookkeeping, via one formula
        (see _formalCharge): charge = (element's neutral valence-electron
        count) − 8 + (final bond-order sum, H included). This single
        formula correctly derives −1 for a carboxylate oxygen, +1 for an
        ammonium nitrogen, and 0 for both pyridine- and pyrrole-type ring
        nitrogens — no element-specific exclusion list needed, because the
        bookkeeping already knows how many bonds (real H included) each
        atom ended up with.

   WHY EXPLICIT HYDROGENS ARE REQUIRED. Without real hydrogen atoms in the
   graph, a terminal heavy atom with one heavy neighbour is structurally
   ambiguous: it could be a saturated −CH3/−NH2/−OH group (no deficit at
   all) or genuinely unsaturated (a real double/triple bond) — pure valence
   counting cannot tell these apart (verified empirically: a heavy-atom-only
   nitrile fragment resolves to the WRONG bond getting the multiple bond).
   With real hydrogens present this ambiguity never arises: a −CH3 carbon
   already has all 4 bonds accounted for (3×H + 1 heavy) and is simply never
   a matching candidate. So perceiveChemistry deliberately returns
   `available: false` for structures with no explicit hydrogens, rather than
   silently falling back to the same kind of geometry-based guessing that
   produced the bug this rewrite fixes. Callers must check `available`
   before using the result, and the UI should disable bead-type prediction
   (and anything else depending on this) when it is false.
   =========================================================================== */

// Neutral free-atom valence electron count (periodic table group), used by
// _formalCharge. Element symbols are upper-cased keys throughout this file.
const _VALENCE_ELECTRONS = {
    H: 1, B: 3, C: 4, N: 5, O: 6, F: 7, SI: 4, P: 5, S: 6, CL: 7, BR: 7, I: 7, SE: 6,
};

// Candidate target valences per element, most-common first. An atom's actual
// target is the smallest option that is still >= its raw bond count, so a
// phosphate P (4-5 bonds) picks 5 instead of the default 3, etc.
const _VALENCE_OPTIONS = {
    H: [1], B: [3, 4], C: [4], N: [3, 4], O: [2, 1, 3],
    F: [1], SI: [4], P: [3, 5], S: [2, 4, 6], CL: [1], BR: [1], I: [1], SE: [2, 4, 6],
};

const _ORGANIC = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'CL', 'BR', 'I']);

/**
 * Formal charge implied by an atom's fully-resolved bond-order sum.
 * charge = (neutral valence electrons) - 8 + (total bond order, H included).
 * H, B, and hypervalent P/S (valence 5 / 6, their normal hypervalent state)
 * get the same special-cased formula xyz2mol uses, since the general octet
 * formula doesn't apply to them.
 * @param {string} el - upper-cased element symbol
 * @param {number} boSum - total resolved bond order (all bonds, H included)
 * @returns {number} formal charge
 */
function _formalCharge(el, boSum) {
    if (el === 'H') return 1 - boSum;
    if (el === 'B') return 3 - boSum;
    if (el === 'P' && boSum === 5) return 0;
    if (el === 'S' && boSum === 6) return 0;
    const ve = _VALENCE_ELECTRONS[el];
    return ve === undefined ? 0 : ve - 8 + boSum;
}

/**
 * Maximum matching on a small general graph, found by backtracking. Returns
 * the largest possible set of edges with no shared endpoint. Small enough
 * graphs (a handful to a few dozen atoms — single organic fragments/rings)
 * make exhaustive backtracking fast; this is not meant for large graphs.
 * @param {Array} nodes - node ids
 * @param {Array<[any,any]>} edgeList - candidate edges between eligible nodes
 * @returns {Array<[any,any]>} the matched pairs
 */
function _maxMatching(nodes, edgeList) {
    const adj = new Map(nodes.map((n) => [n, []]));
    for (const [a, b] of edgeList) { adj.get(a).push(b); adj.get(b).push(a); }

    let best = [];
    function backtrack(remaining, current) {
        if (current.length > best.length) best = current.slice();
        if (remaining.size === 0) return;
        const node = remaining.values().next().value;
        const withoutNode = new Set(remaining);
        withoutNode.delete(node);
        for (const nb of adj.get(node)) {
            if (!withoutNode.has(nb)) continue;
            const next = new Set(withoutNode);
            next.delete(nb);
            backtrack(next, [...current, [node, nb]]);
        }
        backtrack(withoutNode, current); // leave `node` unmatched
    }
    backtrack(new Set(nodes), []);
    return best;
}

/**
 * Core valence-budget + matching engine (see the module-level algorithm
 * comment). Resolves bond orders for one connected component of a molecular
 * graph that includes explicit hydrogen atoms.
 * @param {Map<number,string>} elements - atomIdx -> upper-cased element
 * @param {Map<number,number>} knownCharge - atomIdx -> already-known formal
 *   charge (e.g. from the file), default 0 for anything absent
 * @param {Array<[number,number]>} edges - bonded atom-index pairs
 * @returns {{ order: Map<string,number>, charge: Map<number,number> }}
 *   order keyed "min-max", charge per atom index (all bonds counted, H
 *   included — callers needing "real H count" should count H-element
 *   neighbours directly, it's already explicit in the input graph).
 */
function _resolveBondOrders(elements, knownCharge, edges) {
    const adj = new Map();
    for (const idx of elements.keys()) adj.set(idx, []);
    const order = new Map();
    for (const [a, b] of edges) {
        order.set(`${Math.min(a, b)}-${Math.max(a, b)}`, 1);
        adj.get(a).push(b);
        adj.get(b).push(a);
    }

    const boSum = (idx) => {
        let s = 0;
        for (const nb of adj.get(idx)) s += order.get(`${Math.min(idx, nb)}-${Math.max(idx, nb)}`);
        return s;
    };
    const targetValence = (idx) => {
        const el = elements.get(idx);
        const charge = knownCharge.get(idx) || 0;
        const options = _VALENCE_OPTIONS[el] || [4];
        const rawDegree = adj.get(idx).length;
        const base = options.find((v) => v >= rawDegree) ?? options[options.length - 1];
        return base + charge;
    };
    const deficit = (idx) => targetValence(idx) - boSum(idx);

    // Single matching round: every atom with remaining valence budget pairs
    // with at most one equally-eligible neighbour, becoming a double bond.
    const eligible = [...elements.keys()].filter((idx) => deficit(idx) > 0);
    const eligibleSet = new Set(eligible);
    const pairEdges = edges.filter(([a, b]) => eligibleSet.has(a) && eligibleSet.has(b));
    const matching = _maxMatching(eligible, pairEdges);
    for (const [a, b] of matching) {
        const key = `${Math.min(a, b)}-${Math.max(a, b)}`;
        order.set(key, order.get(key) + 1);
    }

    // Triple-bond escalation: a matched pair that's each still short by
    // exactly 1, with no other eligible neighbour on either side, can only
    // mean a triple bond (nitriles, alkynes) — there's nowhere else for
    // that remaining deficit to go.
    for (const [a, b] of matching) {
        if (deficit(a) === 1 && deficit(b) === 1) {
            const aOthers = adj.get(a).filter((n) => n !== b && deficit(n) > 0);
            const bOthers = adj.get(b).filter((n) => n !== a && deficit(n) > 0);
            if (aOthers.length === 0 && bOthers.length === 0) {
                const key = `${Math.min(a, b)}-${Math.max(a, b)}`;
                order.set(key, order.get(key) + 1);
            }
        }
    }

    const charge = new Map();
    for (const idx of elements.keys()) charge.set(idx, _formalCharge(elements.get(idx), boSum(idx)));

    return { order, charge };
}

/**
 * Find each connected component of a graph (Set of node ids per component).
 * Bond-order resolution runs independently per component, so unrelated
 * molecules sharing one file (e.g. a ligand plus a separately-listed water)
 * never interfere with each other.
 */
function _connectedComponents(nodes, adj) {
    const seen = new Set();
    const components = [];
    for (const start of nodes) {
        if (seen.has(start)) continue;
        const comp = new Set([start]);
        seen.add(start);
        const queue = [start];
        while (queue.length) {
            const u = queue.shift();
            for (const v of (adj.get(u) || [])) {
                if (!seen.has(v)) { seen.add(v); comp.add(v); queue.push(v); }
            }
        }
        components.push(comp);
    }
    return components;
}

/**
 * Perceive whole-molecule chemistry: bond orders, ring/aromaticity, and
 * formal charge, from connectivity + explicit hydrogens alone (see the
 * module-level comment for the algorithm and why explicit H is required).
 *
 * Returns `{ available: false }` when the structure has no explicit
 * hydrogen atoms — callers MUST check `available` before using anything
 * else in the result, and should disable any feature that depends on this
 * (bead-type prediction, whole-molecule SMILES, viewer bond reflection)
 * rather than silently falling back to a guess.
 *
 * @param {object} structure - NGL-style structure (eachAtom/getAtomProxy)
 * @returns {{
 *   available: boolean,
 *   ringAtoms: Set<number>,        // every atom in any ring, aromatic or not
 *   aromaticAtoms: Set<number>,    // ring atoms in a fully-resolved alternating system
 *   branchAtoms: Set<number>,      // heavy atoms with >=3 heavy-atom neighbours
 *   bondOrders: Map<string,number>,// "minIdx-maxIdx" -> order, ALL bonds
 *   charges: Map<number,number>,   // atomIdx -> derived formal charge
 *   hNeighbors: Map<number,number>,// atomIdx -> count of explicit H neighbours
 * }}
 */
export function perceiveChemistry(structure) {
    const empty = {
        available: false, ringAtoms: new Set(), aromaticAtoms: new Set(), branchAtoms: new Set(),
        bondOrders: new Map(), charges: new Map(), hNeighbors: new Map(),
    };
    if (!structure || typeof structure.eachAtom !== 'function') return empty;
    if (!structureHasHydrogens(structure)) return empty;

    const element = new Map();
    structure.eachAtom((a) => element.set(a.index, (a.element || 'C').toUpperCase()));

    // Full graph, hydrogens included as real nodes — see module comment.
    const adj = new Map();
    const edgeSeen = new Set();
    const edges = [];
    const link = (u, v) => { if (!adj.has(u)) adj.set(u, []); adj.get(u).push(v); };
    structure.eachAtom((atom) => {
        const idx = atom.index;
        if (typeof atom.eachBond !== 'function') return;
        atom.eachBond((bond) => {
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const other = i1 === idx ? i2 : (i2 === idx ? i1 : -1);
            if (other < 0) return;
            const a = Math.min(idx, other), b = Math.max(idx, other);
            const key = `${a}-${b}`;
            if (edgeSeen.has(key)) return;
            edgeSeen.add(key);
            edges.push([a, b]);
            link(a, b); link(b, a);
        });
    });

    // Resolve bond orders + charge per connected component independently.
    const bondOrders = new Map();
    const charges = new Map();
    const knownCharge = new Map();
    structure.eachAtom((a) => { if (a.formalCharge) knownCharge.set(a.index, a.formalCharge); });
    for (const comp of _connectedComponents([...element.keys()], adj)) {
        const compElements = new Map([...element].filter(([idx]) => comp.has(idx)));
        const compEdges = edges.filter(([a, b]) => comp.has(a));
        const { order, charge } = _resolveBondOrders(compElements, knownCharge, compEdges);
        for (const [k, v] of order) bondOrders.set(k, v);
        for (const [k, v] of charge) charges.set(k, v);
    }

    const hNeighbors = new Map();
    for (const idx of element.keys()) {
        if (element.get(idx) === 'H') continue;
        let count = 0;
        for (const nb of (adj.get(idx) || [])) if (element.get(nb) === 'H') count += 1;
        hNeighbors.set(idx, count);
    }

    // Ring detection on the heavy-atom subgraph only (BFS spanning tree,
    // handles fused rings correctly — every non-tree edge directly gives a
    // fundamental cycle via its tree-path LCA, unlike DFS back-edges, which
    // miss cycles in polycyclic systems where the closing atom of a second
    // ring is already marked "done" by the time the first ring's DFS
    // reaches it).
    const heavyAdj = new Map();
    for (const idx of element.keys()) {
        if (element.get(idx) === 'H') continue;
        heavyAdj.set(idx, (adj.get(idx) || []).filter((n) => element.get(n) !== 'H'));
    }

    // Branch points: atoms with >=3 heavy-atom neighbours, used by the bead
    // size-class rule (a 4-heavy-atom ring or branched group sizes as S
    // instead of R). Computed from the full heavy-atom adjacency before any
    // ring-pruning below, since branching is a property of the real
    // molecule, not of the pruned ring-only subgraph.
    const branchAtoms = new Set();
    for (const [idx, neighbors] of heavyAdj) {
        if (neighbors.length >= 3) branchAtoms.add(idx);
    }

    const degree = new Map([...heavyAdj].map(([k, v]) => [k, v.length]));
    const removed = new Set();
    let changed = true;
    while (changed) {
        changed = false;
        for (const [k, d] of degree) {
            if (removed.has(k) || d > 1) continue;
            removed.add(k);
            changed = true;
            for (const nb of heavyAdj.get(k)) if (!removed.has(nb)) degree.set(nb, degree.get(nb) - 1);
        }
    }
    const ringAdj = new Map();
    for (const k of heavyAdj.keys()) {
        if (removed.has(k)) continue;
        ringAdj.set(k, heavyAdj.get(k).filter((n) => !removed.has(n)));
    }

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
                for (const v of ringAdj.get(u)) if (!bfsParent.has(v)) { bfsParent.set(v, u); bfsQ.push(v); }
            }
        }
        const treeSet = new Set();
        for (const [n, p] of bfsParent) if (p >= 0) treeSet.add(`${Math.min(n, p)}-${Math.max(n, p)}`);
        const seenEdge = new Set();
        for (const u of ringAdj.keys()) {
            for (const v of ringAdj.get(u)) {
                if (v <= u) continue;
                const ekey = `${u}-${v}`;
                if (treeSet.has(ekey) || seenEdge.has(ekey)) continue;
                seenEdge.add(ekey);
                const pathU = [], pathV = [];
                for (let x = u; x >= 0; x = bfsParent.get(x) ?? -1) pathU.push(x);
                for (let x = v; x >= 0; x = bfsParent.get(x) ?? -1) pathV.push(x);
                const setU = new Map(pathU.map((n, i) => [n, i]));
                let lca = -1, vIdx = -1;
                for (let i = 0; i < pathV.length; i++) if (setU.has(pathV[i])) { lca = pathV[i]; vIdx = i; break; }
                if (lca < 0) continue;
                const cycle = [...pathU.slice(0, setU.get(lca) + 1), ...pathV.slice(0, vIdx).reverse()];
                if (cycle.length >= 3) cycles.push(cycle);
            }
        }
    }

    const ringAtoms = new Set();
    const aromaticAtoms = new Set();
    const isHetero = (idx) => ['N', 'O', 'S', 'SE'].includes(element.get(idx));
    for (const cyc of cycles) {
        if (cyc.length < 4 || cyc.length > 7) continue;
        for (const idx of cyc) ringAtoms.add(idx);

        // Aromatic iff every ring member either has a resolved double bond to
        // ANOTHER ring member, or is a heteroatom with no double bond at all
        // (its lone pair completes the ring instead) — a direct read of the
        // already-resolved bond pattern, not a separate geometric guess.
        const cycSet = new Set(cyc);
        const fullyConjugated = cyc.every((idx) => {
            const hasRingDoubleBond = (heavyAdj.get(idx) || []).some((nb) => {
                if (!cycSet.has(nb)) return false;
                const key = `${Math.min(idx, nb)}-${Math.max(idx, nb)}`;
                return bondOrders.get(key) === 2;
            });
            if (hasRingDoubleBond) return true;
            return isHetero(idx) && !(heavyAdj.get(idx) || []).some((nb) => {
                const key = `${Math.min(idx, nb)}-${Math.max(idx, nb)}`;
                return bondOrders.get(key) === 2;
            });
        });
        if (fullyConjugated) for (const idx of cyc) aromaticAtoms.add(idx);
    }

    return { available: true, ringAtoms, aromaticAtoms, branchAtoms, bondOrders, charges, hNeighbors };
}

/**
 * Build a SMILES string for a set of atoms, given the whole-molecule
 * chemistry already resolved by perceiveChemistry. Used both for one bead's
 * fragment (the common case — bonds leaving the fragment are capped with
 * hydrogen, same as AutoMartini's fragment treatment) and, by passing every
 * heavy atom in the structure, for a whole-molecule SMILES.
 *
 * aromaticNotation=true writes aromatic atoms lowercase and skips Kekulé
 * bond symbols, for matching AutoMartini's open-chain aromatic table keys
 * ("cc", "ccc", ...). Never feed that output to RDKit — it rejects
 * open-chain aromatic SMILES; aromaticNotation=false (the default) produces
 * a normal Kekulé SMILES RDKit can parse.
 *
 * @param {Array} beadAtoms - atom proxies to include (heavy atoms only;
 *   hydrogens are represented implicitly via chemistry.hNeighbors)
 * @param {object} chemistry - perceiveChemistry's result (must have
 *   `available: true`)
 * @param {object} [opts]
 * @param {boolean} [opts.aromaticNotation]
 * @param {number} [opts.startIndex] - atom index to start the SMILES from
 * @returns {string|null}
 */
export function fragmentToSmiles(beadAtoms, chemistry, { aromaticNotation = false, startIndex = null } = {}) {
    const heavy = beadAtoms.filter((a) => (a.element || 'C').toUpperCase() !== 'H');
    if (heavy.length === 0 || !chemistry || !chemistry.available) return null;

    const atomSet = new Set(heavy.map((a) => a.index));
    const { bondOrders, charges, hNeighbors, aromaticAtoms } = chemistry;

    // Collect each internal bond once, reading its already-resolved order.
    const edges = new Map();
    for (const atom of heavy) {
        if (typeof atom.eachBond !== 'function') continue;
        atom.eachBond((bond) => {
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const otherIdx = i1 === atom.index ? i2 : (i2 === atom.index ? i1 : -1);
            if (otherIdx < 0 || !atomSet.has(otherIdx)) return;
            const a = Math.min(atom.index, otherIdx), b = Math.max(atom.index, otherIdx);
            const key = `${a}-${b}`;
            if (edges.has(key)) return;
            const aromatic = aromaticAtoms.has(a) && aromaticAtoms.has(b);
            const order = aromaticNotation && aromatic ? 1 : (bondOrders.get(key) ?? 1);
            edges.set(key, { a, b, order, aromatic });
        });
    }

    const data = new Map();
    for (const atom of heavy) {
        data.set(atom.index, {
            el: (atom.element || 'C').toUpperCase(),
            charge: charges.get(atom.index) ?? 0,
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
    // hCount = real explicit H neighbours, plus one per bond leaving the
    // fragment to ANOTHER HEAVY atom (capped as hydrogen, same
    // fragment-capping convention as before — explicit H neighbours are
    // already counted via hNeighbors and must not be double-counted here).
    for (const atom of heavy) {
        const d = data.get(atom.index);
        const realH = hNeighbors.get(atom.index) ?? 0;
        let externalHeavyBonds = 0;
        if (typeof atom.eachBond === 'function') {
            const structure = atom.structure;
            atom.eachBond((bond) => {
                const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
                const otherIdx = i1 === atom.index ? i2 : (i2 === atom.index ? i1 : -1);
                if (otherIdx < 0 || atomSet.has(otherIdx)) return;
                const other = structure.getAtomProxy(otherIdx);
                if ((other.element || '').toUpperCase() !== 'H') externalHeavyBonds += 1;
            });
        }
        d.hCount = realH + externalHeavyBonds;
    }

    const visited = new Set();
    const inStack = new Set();
    const seenEdges = new Set();
    const ringClosures = new Map();
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
                if (!ringClosures.has(idx)) ringClosures.set(idx, []);
                ringClosures.get(toIdx).push({ digit, writeBond: false });
                ringClosures.get(idx).push({ digit, order, writeBond: true });
            } else if (!visited.has(toIdx)) {
                findBackEdges(toIdx, idx);
            }
        }
        inStack.delete(idx);
    }

    const startIdx = (startIndex != null && atomSet.has(startIndex)) ? startIndex : heavy[0].index;
    findBackEdges(startIdx, -1);
    visited.clear();

    function bondChar(order) {
        return order === 2 ? '=' : order === 3 ? '#' : '';
    }

    function atomToken(idx) {
        const d = data.get(idx);
        if (aromaticNotation && aromaticAtoms.has(idx)) {
            return d.el.charAt(0).toLowerCase() + d.el.slice(1).toLowerCase();
        }
        const sym = d.el === 'CL' ? 'Cl' : d.el === 'BR' ? 'Br' : d.el.charAt(0) + d.el.slice(1).toLowerCase();
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
        const children = (data.get(idx)?.internalBonds ?? []).filter(({ toIdx }) => !visited.has(toIdx));
        if (children.length === 0) return smi;
        for (let i = 0; i < children.length - 1; i++) {
            const { toIdx, order } = children[i];
            if (visited.has(toIdx)) continue;
            smi += `(${bondChar(order)}${dfs(toIdx)})`;
        }
        const last = children[children.length - 1];
        if (!visited.has(last.toIdx)) smi += bondChar(last.order) + dfs(last.toIdx);
        return smi;
    }

    return dfs(startIdx);
}

/**
 * Whole-molecule SMILES, built from every heavy atom in the structure using
 * the same resolved chemistry as the per-bead fragments. Returns null when
 * `chemistry.available` is false (no explicit hydrogens) or the structure
 * has more than one connected component (use per-residue/fragment calls
 * instead for multi-molecule files).
 * @param {object} structure - NGL-style structure
 * @param {object} chemistry - perceiveChemistry's result
 * @returns {string|null}
 */
export function moleculeToSmiles(structure, chemistry) {
    if (!structure || !chemistry || !chemistry.available) return null;

    // Collect indices during the eachAtom pass, then build a fresh proxy per
    // atom — eachAtom's callback argument is a single reused/mutable proxy
    // in real NGL structures, so storing it directly would leave every
    // entry aliased to whichever atom the iteration finished on.
    const heavyIndices = [];
    structure.eachAtom((a) => { if ((a.element || 'C').toUpperCase() !== 'H') heavyIndices.push(a.index); });
    const heavy = heavyIndices.map((idx) => structure.getAtomProxy(idx));

    // A multi-molecule file (e.g. a ligand plus a separately-listed ion)
    // would otherwise silently produce a SMILES for just whichever fragment
    // the DFS happens to start in, rather than the null this is documented
    // to return.
    const heavySet = new Set(heavyIndices);
    const adj = new Map(heavyIndices.map((idx) => [idx, []]));
    for (const atom of heavy) {
        if (typeof atom.eachBond !== 'function') continue;
        atom.eachBond((bond) => {
            const other = bond.atomIndex1 === atom.index ? bond.atomIndex2 : bond.atomIndex1;
            if (heavySet.has(other)) adj.get(atom.index).push(other);
        });
    }
    if (_connectedComponents(heavyIndices, adj).length > 1) return null;

    return fragmentToSmiles(heavy, chemistry);
}

/**
 * Number of H-bond donor atoms (N/O/S with at least one real bonded
 * hydrogen) in a bead, read directly from perceiveChemistry's hNeighbors —
 * always the real, explicit count, since prediction is gated on explicit
 * hydrogens being present in the first place.
 * @param {object} bead - has an `atoms` array of atom proxies
 * @param {object} chemistry - perceiveChemistry's result
 * @returns {number}
 */
export function beadDonorCount(bead, chemistry) {
    let donors = 0;
    for (const atom of bead.atoms) {
        const el = (atom.element || '').toUpperCase();
        if ((el === 'N' || el === 'O' || el === 'S') && (chemistry.hNeighbors.get(atom.index) ?? 0) > 0) {
            donors += 1;
        }
    }
    return donors;
}

// Heteroatoms whose H-capped representation (see fragmentToSmiles) can read
// as a different, more terminal group than the real structure — an ether
// capped at the O looks like an alcohol, a secondary amine capped at the N
// looks like primary, a thioether like a thiol. Carbon doesn't have this
// problem (a capped alkyl chain still reads as an alkyl chain), so it's
// excluded.
const _CAPPING_SENSITIVE = new Set(['N', 'O', 'S', 'P']);

/**
 * Heteroatoms (N/O/S/P) in a bead that have at least one bond leaving the
 * bead to another heavy atom — i.e. a bead boundary that cuts through a
 * functional group rather than containing it whole. Purely structural (only
 * needs atom.eachBond), independent of perceiveChemistry/chemistry.available,
 * so it stays useful even without explicit hydrogens. Used to surface a
 * "this bead may be chopping a chemical group" warning.
 * @param {Array} beadAtoms - atom proxies in the bead
 * @returns {Array} the offending atom proxies (empty if none)
 */
export function cappedHeteroatoms(beadAtoms) {
    const atomSet = new Set(beadAtoms.map((a) => a.index));
    const result = [];
    for (const atom of beadAtoms) {
        const el = (atom.element || '').toUpperCase();
        if (!_CAPPING_SENSITIVE.has(el) || typeof atom.eachBond !== 'function') continue;
        const structure = atom.structure;
        let capped = false;
        atom.eachBond((bond) => {
            const i1 = bond.atomIndex1, i2 = bond.atomIndex2;
            const otherIdx = i1 === atom.index ? i2 : (i2 === atom.index ? i1 : -1);
            if (otherIdx < 0 || atomSet.has(otherIdx)) return;
            const other = structure.getAtomProxy(otherIdx);
            if ((other.element || '').toUpperCase() !== 'H') capped = true;
        });
        if (capped) result.push(atom);
    }
    return result;
}

/**
 * Count distinct residues in a structure (by resno+chain, falling back to
 * resno alone), for surfacing a "this structure has more than one residue"
 * warning. Loaded files are expected to be a single molecule for mapping
 * purposes; multiple residues usually means solvent, ions, or multiple
 * copies came along with the structure.
 * @param {object} structure - NGL-style structure
 * @returns {number}
 */
export function countResidues(structure) {
    const seen = new Set();
    if (structure && typeof structure.eachAtom === 'function') {
        structure.eachAtom((atom) => {
            seen.add(`${atom.chainname ?? ''}/${atom.resno}/${atom.resname ?? ''}`);
        });
    }
    return seen.size;
}

// Elements in a period higher than the 3rd (Br/Se = period 4, I = period 5).
// Per the Martini 3 SI's default bead-size convention, these count as TWO
// non-hydrogen atoms when sizing a bead, since they're physically bulkier
// than a typical 2nd/3rd-period heavy atom (S/P/Cl/Si stay at normal weight
// — the SI's example is iodine). Used both for whole-molecule heavy-atom
// counts (the Mismatch panel) and per-bead weighted counts (bead size-class
// prediction), so both stay on the same scale.
const _PERIOD_4_PLUS = new Set(['BR', 'SE', 'I']);

/**
 * Bead-sizing weight for one element: 2 for a period->=4 atom, 1 otherwise.
 * @param {string} element - element symbol, any case
 * @returns {number} 1 or 2
 */
export function heavyAtomWeight(element) {
    return _PERIOD_4_PLUS.has((element || '').toUpperCase()) ? 2 : 1;
}

/**
 * Sum of heavyAtomWeight() over every non-hydrogen atom in a structure —
 * the weighted heavy-atom total used as the Mismatch panel's reference
 * value, kept on the same scale as bead-type-implied expected counts.
 * @param {object} structure - NGL-style structure (eachAtom)
 * @returns {number}
 */
export function weightedHeavyAtomCount(structure) {
    let total = 0;
    if (structure && typeof structure.eachAtom === 'function') {
        structure.eachAtom((atom) => {
            const el = (atom.element || '').toUpperCase();
            if (el !== 'H') total += heavyAtomWeight(el);
        });
    }
    return total;
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
