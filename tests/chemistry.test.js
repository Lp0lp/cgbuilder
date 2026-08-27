import { describe, it, expect } from 'vitest';
import {
    perceiveChemistry, fragmentToSmiles, moleculeToSmiles, beadDonorCount, countResidues,
    heavyAtomWeight, weightedHeavyAtomCount,
} from '../scripts/chemistry.js';
import { buildStructure } from './helpers/mockStructure.js';

// Builds an n-membered ring of `ringElement` atoms, each carrying `hPerAtom`
// explicit hydrogens except where overridden in `noH` (a set of ring
// positions, 0-indexed) -- the common shape for the aromatic/ring fixtures
// below. Ring atom indices are 0..n-1; any hydrogens come after.
function ringWithH(ringElements, hPerAtom, noH = new Set()) {
    const n = ringElements.length;
    const atomDefs = ringElements.map((element) => ({ element }));
    const bonds = [];
    for (let i = 0; i < n; i++) bonds.push({ a: i, b: (i + 1) % n });
    for (let i = 0; i < n; i++) {
        if (noH.has(i)) continue;
        for (let h = 0; h < hPerAtom; h++) {
            atomDefs.push({ element: 'H' });
            bonds.push({ a: i, b: atomDefs.length - 1 });
        }
    }
    return buildStructure(atomDefs, bonds);
}

function heavyAtomsOf(structure) {
    const out = [];
    structure.eachAtom((a) => { if (a.element !== 'H') out.push(a); });
    return out;
}

describe('perceiveChemistry — gating on explicit hydrogens', () => {
    it('is unavailable for a structure with no explicit hydrogens', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'O' }], [{ a: 0, b: 1 }]);
        const chem = perceiveChemistry(structure);
        expect(chem.available).toBe(false);
    });

    it('is available once at least one explicit hydrogen is present', () => {
        const structure = buildStructure(
            [{ element: 'C' }, { element: 'H' }], [{ a: 0, b: 1 }],
        );
        expect(perceiveChemistry(structure).available).toBe(true);
    });

    it('returns the empty/unavailable shape for a missing or malformed structure', () => {
        expect(perceiveChemistry(null).available).toBe(false);
        expect(perceiveChemistry({}).available).toBe(false);
    });
});

describe('perceiveChemistry — ring detection and aromaticity', () => {
    it('flags benzene (6 CH) as fully aromatic', () => {
        const structure = ringWithH(['C', 'C', 'C', 'C', 'C', 'C'], 1);
        const chem = perceiveChemistry(structure);
        expect(chem.ringAtoms.size).toBe(6);
        expect(chem.aromaticAtoms.size).toBe(6);
    });

    it('excludes the furan oxygen from any double bond, but still counts it aromatic', () => {
        // Ring: O(0)-C(1)-C(2)-C(3)-C(4), O has no H, each carbon has 1.
        const structure = ringWithH(['O', 'C', 'C', 'C', 'C'], 1, new Set([0]));
        const chem = perceiveChemistry(structure);
        expect(chem.aromaticAtoms.size).toBe(5);
        expect(chem.bondOrders.get('0-1')).toBe(1);
        expect(chem.bondOrders.get('0-4')).toBe(1);
        // exactly two ring C=C double bonds among the four carbons
        const ringDoubles = ['1-2', '2-3', '3-4'].filter((k) => chem.bondOrders.get(k) === 2);
        expect(ringDoubles).toHaveLength(2);
    });

    it('excludes the thiophene sulfur from any double bond the same way', () => {
        const structure = ringWithH(['S', 'C', 'C', 'C', 'C'], 1, new Set([0]));
        const chem = perceiveChemistry(structure);
        expect(chem.aromaticAtoms.size).toBe(5);
        expect(chem.bondOrders.get('0-1')).toBe(1);
        expect(chem.bondOrders.get('0-4')).toBe(1);
    });

    it('gives pyridine-type N (no H) a double bond, unlike furan/thiophene', () => {
        const structure = ringWithH(['N', 'C', 'C', 'C', 'C', 'C'], 1, new Set([0]));
        const chem = perceiveChemistry(structure);
        expect(chem.aromaticAtoms.size).toBe(6);
        const nHasDouble = chem.bondOrders.get('0-1') === 2 || chem.bondOrders.get('0-5') === 2;
        expect(nHasDouble).toBe(true);
        expect(chem.charges.get(0)).toBe(0);
    });

    it('gives pyrrole-type N (with H) zero double bonds and zero charge', () => {
        const structure = ringWithH(['N', 'C', 'C', 'C', 'C'], 1); // N keeps its H like every other ring atom
        const chem = perceiveChemistry(structure);
        expect(chem.aromaticAtoms.size).toBe(5);
        expect(chem.bondOrders.get('0-1')).toBe(1);
        expect(chem.bondOrders.get('0-4')).toBe(1);
        expect(chem.charges.get(0)).toBe(0);
    });

    it('does not flag a saturated ring (cyclohexane) as aromatic, but still counts it as a ring', () => {
        const structure = ringWithH(['C', 'C', 'C', 'C', 'C', 'C'], 2); // 2 H per carbon -> fully saturated
        const chem = perceiveChemistry(structure);
        expect(chem.ringAtoms.size).toBe(6);
        expect(chem.aromaticAtoms.size).toBe(0);
        for (let i = 0; i < 6; i++) {
            const next = (i + 1) % 6;
            expect(chem.bondOrders.get(`${Math.min(i, next)}-${Math.max(i, next)}`)).toBe(1);
        }
    });

    it('finds aromaticity in both rings of a fused bicyclic system (benzofuran — regression for the DFS->BFS cycle-finder fix)', () => {
        // Benzofuran, standard numbering: O1(0),C2(1),C3(2),C3a(3),C7a(4) form
        // the furan ring; C3a(3)/C7a(4) are the fusion atoms (3 ring bonds,
        // no H); C4(5),C5(6),C6(7),C7(8) are the extra benzo-ring carbons.
        // This is exactly the class of molecule (a fused heteroaromatic ring)
        // that motivated the original DFS->BFS cycle-finder fix.
        const atomDefs = [
            { element: 'O' }, { element: 'C' }, { element: 'C' }, { element: 'C' }, { element: 'C' },
            { element: 'C' }, { element: 'C' }, { element: 'C' }, { element: 'C' },
        ];
        const bonds = [
            { a: 0, b: 1 }, { a: 1, b: 2 }, { a: 2, b: 3 }, { a: 3, b: 4 }, { a: 4, b: 0 }, // furan ring
            { a: 3, b: 5 }, { a: 5, b: 6 }, { a: 6, b: 7 }, { a: 7, b: 8 }, { a: 8, b: 4 }, // benzo ring
        ];
        for (const idx of [1, 2, 5, 6, 7, 8]) { // CH positions; O(0) and the two fusion carbons get no H
            atomDefs.push({ element: 'H' });
            bonds.push({ a: idx, b: atomDefs.length - 1 });
        }
        const structure = buildStructure(atomDefs, bonds);
        const chem = perceiveChemistry(structure);
        expect(chem.ringAtoms.size).toBe(9);
        expect(chem.aromaticAtoms.size).toBe(9);
        expect(chem.charges.get(0)).toBe(0);

        // Verified separately via RDKit to canonicalize identically to a
        // textbook benzofuran SMILES ("c1ccc2occc2c1").
        const heavy = heavyAtomsOf(structure);
        expect(fragmentToSmiles(heavy, chem)).toBe('O1(C=CC2(=C1C=CC=C2))');
    });
});

describe('perceiveChemistry — bond order and formal charge resolution', () => {
    it('resolves an amide: carbonyl C=O wins, C-N stays single, N gets its 2 H', () => {
        // C(0)=O(1), C(0)-N(2), N(2)-H(3), N(2)-H(4)
        const structure = buildStructure(
            [{ element: 'C' }, { element: 'O' }, { element: 'N' }, { element: 'H' }, { element: 'H' }],
            [{ a: 0, b: 1 }, { a: 0, b: 2 }, { a: 2, b: 3 }, { a: 2, b: 4 }],
        );
        const chem = perceiveChemistry(structure);
        expect(chem.bondOrders.get('0-1')).toBe(2);
        expect(chem.bondOrders.get('0-2')).toBe(1);
        expect(chem.charges.get(2)).toBe(0);
    });

    it('derives -1 on the non-double-bonded oxygen of a deprotonated carboxylate', () => {
        // H3C(3)-C(4)(=O(5))-O(6), no H on either oxygen.
        const structure = buildStructure(
            [{ element: 'H' }, { element: 'H' }, { element: 'H' }, { element: 'C' }, { element: 'C' }, { element: 'O' }, { element: 'O' }],
            [{ a: 0, b: 3 }, { a: 1, b: 3 }, { a: 2, b: 3 }, { a: 3, b: 4 }, { a: 4, b: 5 }, { a: 4, b: 6 }],
        );
        const chem = perceiveChemistry(structure);
        const oCharges = [chem.charges.get(5), chem.charges.get(6)].sort();
        expect(oCharges).toEqual([-1, 0]);
        const heavy = heavyAtomsOf(structure);
        expect(fragmentToSmiles(heavy, chem)).toBe('CC(=O)[O-]');
    });

    it('derives +1 on a methylammonium nitrogen', () => {
        const structure = buildStructure(
            [{ element: 'H' }, { element: 'H' }, { element: 'H' }, { element: 'C' }, { element: 'N' }, { element: 'H' }, { element: 'H' }, { element: 'H' }],
            [{ a: 0, b: 3 }, { a: 1, b: 3 }, { a: 2, b: 3 }, { a: 3, b: 4 }, { a: 4, b: 5 }, { a: 4, b: 6 }, { a: 4, b: 7 }],
        );
        const chem = perceiveChemistry(structure);
        expect(chem.charges.get(4)).toBe(1);
        expect(fragmentToSmiles(heavyAtomsOf(structure), chem)).toBe('C[NH3+]');
    });

    it('resolves acetonitrile with a real triple bond on C#N, not on the methyl C-C bond', () => {
        // This is the exact case that fails without explicit hydrogens: pure
        // heavy-atom valence counting can't tell a saturated terminal CH3
        // apart from an unsaturated terminal atom.
        const structure = buildStructure(
            [{ element: 'H' }, { element: 'H' }, { element: 'H' }, { element: 'C' }, { element: 'C' }, { element: 'N' }],
            [{ a: 0, b: 3 }, { a: 1, b: 3 }, { a: 2, b: 3 }, { a: 3, b: 4 }, { a: 4, b: 5 }],
        );
        const chem = perceiveChemistry(structure);
        expect(chem.bondOrders.get('3-4')).toBe(1);
        expect(chem.bondOrders.get('4-5')).toBe(3);
        expect(fragmentToSmiles(heavyAtomsOf(structure), chem)).toBe('CC#N');
    });
});

describe('fragmentToSmiles', () => {
    it('returns null when chemistry is unavailable', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'O' }], [{ a: 0, b: 1 }]);
        const chem = perceiveChemistry(structure); // no explicit H -> unavailable
        expect(fragmentToSmiles(heavyAtomsOf(structure), chem)).toBeNull();
    });

    it('does not bracket a neutral organic atom just because it has implicit H', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'H' }, { element: 'H' }, { element: 'H' }, { element: 'H' }],
            [{ a: 0, b: 1 }, { a: 0, b: 2 }, { a: 0, b: 3 }, { a: 0, b: 4 }]);
        const chem = perceiveChemistry(structure);
        expect(fragmentToSmiles(heavyAtomsOf(structure), chem)).toBe('C');
    });

    it('closes a ring with a digit for a 3-membered carbocycle', () => {
        const atomDefs = [{ element: 'C' }, { element: 'C' }, { element: 'C' }];
        const bonds = [{ a: 0, b: 1 }, { a: 1, b: 2 }, { a: 2, b: 0 }];
        for (let i = 0; i < 3; i++) {
            for (let h = 0; h < 2; h++) {
                atomDefs.push({ element: 'H' });
                bonds.push({ a: i, b: atomDefs.length - 1 });
            }
        }
        const structure = buildStructure(atomDefs, bonds);
        const chem = perceiveChemistry(structure);
        expect(fragmentToSmiles(heavyAtomsOf(structure), chem)).toBe('C1(CC1)');
    });

    it('writes aromatic-notation lowercase atoms when requested', () => {
        const structure = ringWithH(['C', 'C', 'C', 'C', 'C', 'C'], 1);
        const chem = perceiveChemistry(structure);
        const smi = fragmentToSmiles(heavyAtomsOf(structure), chem, { aromaticNotation: true });
        expect(smi).toBe('c1(ccccc1)');
    });
});

describe('moleculeToSmiles', () => {
    it('matches the fragment SMILES for a whole single-molecule structure', () => {
        const structure = ringWithH(['O', 'C', 'C', 'C', 'C'], 1, new Set([0]));
        const chem = perceiveChemistry(structure);
        expect(moleculeToSmiles(structure, chem)).toBe(fragmentToSmiles(heavyAtomsOf(structure), chem));
    });

    it('returns null when chemistry is unavailable', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'O' }], [{ a: 0, b: 1 }]);
        expect(moleculeToSmiles(structure, perceiveChemistry(structure))).toBeNull();
    });
});

describe('beadDonorCount', () => {
    it('counts only N/O/S atoms with at least one real explicit hydrogen', () => {
        // O(0)-H(2): donor. N(1) with no H: not a donor.
        const structure = buildStructure(
            [{ element: 'O' }, { element: 'N' }, { element: 'H' }, { element: 'C' }],
            [{ a: 0, b: 2 }, { a: 0, b: 3 }, { a: 1, b: 3 }],
        );
        const chem = perceiveChemistry(structure);
        const bead = { atoms: [structure.getAtomProxy(0), structure.getAtomProxy(1)] };
        expect(beadDonorCount(bead, chem)).toBe(1);
    });
});

describe('countResidues', () => {
    it('counts a single residue as 1', () => {
        const structure = buildStructure([{ element: 'C', resno: 1 }, { element: 'O', resno: 1 }], [{ a: 0, b: 1 }]);
        expect(countResidues(structure)).toBe(1);
    });

    it('counts multiple distinct residues', () => {
        const structure = buildStructure(
            [{ element: 'C', resno: 1 }, { element: 'O', resno: 2 }, { element: 'N', resno: 3 }], [],
        );
        expect(countResidues(structure)).toBe(3);
    });
});

describe('heavyAtomWeight / weightedHeavyAtomCount', () => {
    it('weights F as 0.5, period-4+ elements as 2, and everything else as 1', () => {
        expect(heavyAtomWeight('C')).toBe(1);
        expect(heavyAtomWeight('F')).toBe(0.5);
        expect(heavyAtomWeight('S')).toBe(1);  // period 3 -- explicitly not weighted
        expect(heavyAtomWeight('Cl')).toBe(1); // period 3
        expect(heavyAtomWeight('Br')).toBe(2);
        expect(heavyAtomWeight('Se')).toBe(2);
        expect(heavyAtomWeight('I')).toBe(2);
    });

    it('sums per-atom weights across a whole structure, excluding hydrogens', () => {
        // 1 C + 1 Br + 1 F (+ hydrogens, which never count): weighted = 1 + 2 + 0.5 = 3.5
        const structure = buildStructure(
            [{ element: 'C' }, { element: 'Br' }, { element: 'F' }, { element: 'H' }, { element: 'H' }, { element: 'H' }],
            [{ a: 0, b: 1 }, { a: 0, b: 2 }, { a: 0, b: 3 }, { a: 0, b: 4 }, { a: 0, b: 5 }],
        );
        expect(weightedHeavyAtomCount(structure)).toBe(3.5);
    });
});

describe('perceiveChemistry — branchAtoms', () => {
    it('flags a branch point (>=3 heavy neighbours) and nothing else', () => {
        // Isobutane-like: central C0 bonded to 3 methyl carbons + 1 H.
        const structure = buildStructure(
            [
                { element: 'C' }, { element: 'C' }, { element: 'C' }, { element: 'C' }, // C0..C3
                { element: 'H' }, // on C0
                { element: 'H' }, { element: 'H' }, { element: 'H' }, // on C1
                { element: 'H' }, { element: 'H' }, { element: 'H' }, // on C2
                { element: 'H' }, { element: 'H' }, { element: 'H' }, // on C3
            ],
            [
                { a: 0, b: 1 }, { a: 0, b: 2 }, { a: 0, b: 3 }, { a: 0, b: 4 },
                { a: 1, b: 5 }, { a: 1, b: 6 }, { a: 1, b: 7 },
                { a: 2, b: 8 }, { a: 2, b: 9 }, { a: 2, b: 10 },
                { a: 3, b: 11 }, { a: 3, b: 12 }, { a: 3, b: 13 },
            ],
        );
        const chem = perceiveChemistry(structure);
        expect(chem.branchAtoms.has(0)).toBe(true);
        expect(chem.branchAtoms.has(1)).toBe(false);
        expect(chem.branchAtoms.has(2)).toBe(false);
        expect(chem.branchAtoms.has(3)).toBe(false);
    });
});
