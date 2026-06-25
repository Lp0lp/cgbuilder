import { describe, it, expect } from 'vitest';
import {
    perceiveChemistry, fragmentToSmiles, realHydrogenCount,
    beadDonorCount, structureHasHydrogens,
} from '../scripts/chemistry.js';
import { buildStructure, regularPolygon, ringBonds } from './helpers/mockStructure.js';

describe('perceiveChemistry', () => {
    it('flags a simple benzene-like ring as aromatic', () => {
        const atomDefs = regularPolygon(6, 1.39, 'C');
        const bonds = ringBonds([0, 1, 2, 3, 4, 5]);
        const structure = buildStructure(atomDefs, bonds);

        const { aromaticAtoms } = perceiveChemistry(structure);
        expect(aromaticAtoms.size).toBe(6);
        for (let i = 0; i < 6; i++) expect(aromaticAtoms.has(i)).toBe(true);
    });

    it('finds aromaticity in both rings of a fused bicyclic system (regression for DFS->BFS fix)', () => {
        // Synthetic 9-atom, 10-bond fused hexagon+pentagon (all-carbon, every
        // ring bond ~1.39 A). Atoms 0,1 are the shared edge. This topology is
        // exactly the shape that broke under the old DFS back-edge cycle
        // finder: the closing atom of the second ring was already "done" by
        // the time the first ring's DFS reached it, so it never got flagged.
        const L = 1.39;
        const s1 = { element: 'C', x: -L / 2, y: 0, z: 0 };
        const s2 = { element: 'C', x: L / 2, y: 0, z: 0 };

        // Hexagon: shares edge s1(idx0)-s2(idx1), extends to +y.
        const hexCenter = { x: 0, y: L * Math.sqrt(3) / 2, z: 0 };
        const hex = [];
        for (let k = 0; k < 4; k++) {
            const theta = (Math.PI / 180) * (60 * k);
            hex.push({
                element: 'C',
                x: hexCenter.x + L * Math.cos(theta),
                y: hexCenter.y + L * Math.sin(theta),
                z: 0,
            });
        }
        // Pentagon: shares the same edge, extends to -y.
        const R5 = L / (2 * Math.sin(Math.PI / 5));
        const penCenter = { x: 0, y: -Math.sqrt(R5 * R5 - (L / 2) * (L / 2)), z: 0 };
        const pen = [];
        for (const deg of [198, 270, 342]) {
            const theta = (Math.PI / 180) * deg;
            pen.push({
                element: 'C',
                x: penCenter.x + R5 * Math.cos(theta),
                y: penCenter.y + R5 * Math.sin(theta),
                z: 0,
            });
        }

        // Indices: 0=s1, 1=s2, 2..5=hex extras, 6..8=pentagon extras.
        const atomDefs = [s1, s2, ...hex, ...pen];
        const bonds = [
            ...ringBonds([0, 1, 2, 3, 4, 5]), // hexagon (s1-s2-V0-V60-V120-V180)
            { a: 0, b: 6 }, { a: 6, b: 7 }, { a: 7, b: 8 }, { a: 8, b: 1 }, // pentagon remainder
        ];
        const structure = buildStructure(atomDefs, bonds);

        const { aromaticAtoms } = perceiveChemistry(structure);
        expect(aromaticAtoms.size).toBe(9);
        for (let i = 0; i < 9; i++) expect(aromaticAtoms.has(i)).toBe(true);
    });

    it('does not flag a plain alkane chain as aromatic or multiply-bonded', () => {
        const atomDefs = [
            { element: 'C', x: 0, y: 0, z: 0 },
            { element: 'C', x: 1.54, y: 0, z: 0 },
            { element: 'C', x: 3.08, y: 0, z: 0 },
        ];
        const bonds = [{ a: 0, b: 1 }, { a: 1, b: 2 }];
        const structure = buildStructure(atomDefs, bonds);

        const { aromaticAtoms, bondOrders } = perceiveChemistry(structure);
        expect(aromaticAtoms.size).toBe(0);
        expect(bondOrders.size).toBe(0);
    });

    it('detects an isolated C=C double bond from geometry', () => {
        const atomDefs = [
            { element: 'C', x: 0, y: 0, z: 0 },
            { element: 'C', x: 1.34, y: 0, z: 0 },
        ];
        const structure = buildStructure(atomDefs, [{ a: 0, b: 1 }]);

        const { bondOrders } = perceiveChemistry(structure);
        expect(bondOrders.get('0-1')).toBe(2);
    });

    it('resolves an amide: C=O wins the double bond, C-N stays single', () => {
        // Real amide geometry: C=O ~1.23 A, C-N ~1.33 A. The naive nearest-
        // reference-length check classifies BOTH as "double" in isolation
        // (1.33 is closer to the C-N double reference of 1.29 than to single
        // 1.47) -- the greedy "most-shortened-first, one double bond per
        // atom" resolution is what's supposed to keep only the carbonyl.
        const atomDefs = [
            { element: 'C', x: 0, y: 0, z: 0 },      // 0: carbonyl C
            { element: 'O', x: 1.23, y: 0, z: 0 },   // 1: carbonyl O
            { element: 'N', x: -1.33, y: 0, z: 0 },  // 2: amide N
        ];
        const structure = buildStructure(atomDefs, [{ a: 0, b: 1 }, { a: 0, b: 2 }]);

        const { bondOrders } = perceiveChemistry(structure);
        expect(bondOrders.get('0-1')).toBe(2);
        expect(bondOrders.has('0-2')).toBe(false); // left single (absent = order 1)
    });
});

describe('fragmentToSmiles', () => {
    it('does not add brackets to a neutral organic atom just because hCount > 0', () => {
        const structure = buildStructure([{ element: 'C' }]);
        const atoms = [];
        structure.eachAtom((a) => atoms.push(a));
        expect(fragmentToSmiles(atoms)).toBe('C');
    });

    it('brackets a charged atom and shows the charge', () => {
        const structure = buildStructure([{ element: 'N', formalCharge: 1 }]);
        const atoms = [];
        structure.eachAtom((a) => atoms.push(a));
        expect(fragmentToSmiles(atoms)).toBe('[NH4+]');
    });

    it('closes a ring with a digit for a 3-membered carbocycle', () => {
        const atomDefs = regularPolygon(3, 1.54, 'C');
        const structure = buildStructure(atomDefs, ringBonds([0, 1, 2]));
        const atoms = [];
        structure.eachAtom((a) => atoms.push(a));
        const smi = fragmentToSmiles(atoms);
        expect(smi).toBe('C1(CC1)');
    });

    it('writes aromatic-notation lowercase atoms when aromaticAtoms is supplied', () => {
        const atomDefs = regularPolygon(6, 1.39, 'C');
        const structure = buildStructure(atomDefs, ringBonds([0, 1, 2, 3, 4, 5]));
        const atoms = [];
        structure.eachAtom((a) => atoms.push(a));
        const aromaticAtoms = new Set([0, 1, 2, 3, 4, 5]);
        const smi = fragmentToSmiles(atoms, { aromaticNotation: true, aromaticAtoms });
        expect(smi).toBe('c1(ccccc1)');
    });
});

describe('realHydrogenCount / beadDonorCount (explicit-H aware)', () => {
    it('trusts the explicit H count when the molecule carries real hydrogens', () => {
        // Amine N with one explicit H neighbour plus two heavy neighbours.
        const atomDefs = [
            { element: 'N' }, { element: 'C' }, { element: 'C' }, { element: 'H' },
        ];
        const bonds = [{ a: 0, b: 1 }, { a: 0, b: 2 }, { a: 0, b: 3 }];
        const structure = buildStructure(atomDefs, bonds);
        const n = structure.getAtomProxy(0);
        expect(realHydrogenCount(n, true)).toBe(1);
    });

    it('reports zero H on a pyridine-type aromatic N even when explicit (ground truth, not estimated)', () => {
        const atomDefs = [
            { element: 'N', aromatic: true }, { element: 'C', aromatic: true }, { element: 'C', aromatic: true },
        ];
        const bonds = [{ a: 0, b: 1, order: 1.5 }, { a: 0, b: 2, order: 1.5 }];
        const structure = buildStructure(atomDefs, bonds);
        const n = structure.getAtomProxy(0);
        expect(realHydrogenCount(n, true)).toBe(0);
    });

    it('estimates pyridine-type N implicit H as 0 when no explicit hydrogens exist at all', () => {
        const atomDefs = [
            { element: 'N', aromatic: true }, { element: 'C', aromatic: true }, { element: 'C', aromatic: true },
        ];
        const bonds = [{ a: 0, b: 1, order: 1.5 }, { a: 0, b: 2, order: 1.5 }];
        const structure = buildStructure(atomDefs, bonds);
        const n = structure.getAtomProxy(0);
        expect(realHydrogenCount(n, false)).toBe(0);
    });

    it('counts donors across a bead using explicit hydrogens', () => {
        const atomDefs = [
            { element: 'O' }, { element: 'C' }, { element: 'H' }, // hydroxyl O-H
            { element: 'N' }, { element: 'C' },                   // tertiary-like N, no H
        ];
        const bonds = [{ a: 0, b: 1 }, { a: 0, b: 2 }, { a: 3, b: 4 }];
        const structure = buildStructure(atomDefs, bonds);
        const bead = { atoms: [structure.getAtomProxy(0), structure.getAtomProxy(3)] };
        expect(beadDonorCount(bead, true)).toBe(1);
    });
});

describe('structureHasHydrogens', () => {
    it('detects explicit hydrogens when present', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'H' }]);
        expect(structureHasHydrogens(structure)).toBe(true);
    });

    it('returns false for a heavy-atom-only structure', () => {
        const structure = buildStructure([{ element: 'C' }, { element: 'O' }]);
        expect(structureHasHydrogens(structure)).toBe(false);
    });
});
