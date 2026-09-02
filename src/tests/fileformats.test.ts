import { describe, it, expect } from 'vitest';
import {
    generateNDX, generateMap, generatePythonAssignments, generateGRO,
    generateBartender, generatePyCGTOOL, parseShakerMapping,
} from '../fileformats.js';
import { EXAMPLE_MAPPING } from '../example.js';
import type { BeadCollection } from '../bead.js';

// Minimal bead/collection mocks -- just the shape each generator reads, cast
// to the nominal BeadCollection at the boundary (the generators only touch a
// handful of fields).
function makeCollection(beads: unknown[], atomNameMap: Record<number, string> = {}): BeadCollection {
    return {
        beads,
        get assignedBeads() { return (beads as { atoms: unknown[] }[]).filter((b) => b.atoms.length > 0); },
        atomName(atom: { index: number }) { return atomNameMap[atom.index] ?? `A${atom.index}`; },
        expandedAtomNames(bead: { atoms: { index: number }[]; expandedAtoms?: () => { index: number }[] }) {
            const atoms = bead.expandedAtoms ? bead.expandedAtoms() : bead.atoms;
            return atoms.map((a) => atomNameMap[a.index] ?? `A${a.index}`);
        },
    } as unknown as BeadCollection;
}

describe('generateNDX', () => {
    it('writes one [ name ] group per bead with 1-indexed atom numbers', () => {
        const collection = makeCollection([
            { name: 'B0', atoms: [{ index: 0 }, { index: 1 }] },
            { name: 'B1', atoms: [{ index: 2 }] },
        ]);
        const ndx = generateNDX(collection);
        expect(ndx).toBe('[ B0 ]\n1 2 \n\n[ B1 ]\n3 \n\n');
    });
});

describe('generateMap', () => {
    it('lists each atom once with the bead names that reference it', () => {
        const collection = makeCollection(
            [{ name: 'BB', atoms: [{ index: 0 }, { index: 1 }] }],
            { 0: 'CA', 1: 'CB' },
        );
        const map = generateMap(collection);
        expect(map).toContain('[ to ]\nmartini');
        expect(map).toContain('BB ');
        expect(map).toContain('[ atoms ]');
        expect(map).toMatch(/1\tCA\tBB/);
        expect(map).toMatch(/2\tCB\tBB/);
    });
});

describe('generatePythonAssignments', () => {
    it('emits a Shaker-style mapping dict literal', () => {
        const collection = makeCollection(
            [{ name: 'BB', type: 'P2', charge: 0, resname: 'UNK', atoms: [{ index: 0 }] }],
            { 0: 'CA' },
        );
        const out = generatePythonAssignments(collection);
        expect(out).toContain('mapping = {');
        expect(out).toContain('"UNK": {');
        expect(out).toContain('"BB": {"type": "P2", "charge": 0, "atoms": [\'CA\']}');
    });

    it('returns an empty string when there are no beads', () => {
        expect(generatePythonAssignments(makeCollection([]))).toBe('');
    });
});

describe('generateGRO', () => {
    it('converts bead centres from Angstrom to nm and pads fixed-width columns', () => {
        const collection = makeCollection([
            { name: 'BB', resname: 'UNK', resid: 1, center: { x: 10, y: 20, z: 30 }, atoms: [{}] },
        ]);
        const gro = generateGRO(collection);
        const lines = gro.split('\n');
        expect(lines[0]).toBe('Generated with cgbuilder');
        expect(lines[1]).toBe('1');
        expect(lines[2]).toContain('1.000');  // 10 A -> 1.000 nm
        expect(lines[2]).toContain('2.000');
        expect(lines[2]).toContain('3.000');
        expect(gro).toContain('10 10 10');
    });
});

describe('generateBartender', () => {
    it('writes BEADS header and one line per bead with 1-based atom indices', () => {
        const b0 = [{ index: 0 }, { index: 1 }];
        const b1 = [{ index: 2 }];
        const collection = makeCollection([
            { name: 'B0', atoms: b0, expandedAtoms: () => b0 },
            { name: 'B1', atoms: b1, expandedAtoms: () => b1 },
        ]);
        expect(generateBartender(collection)).toBe('BEADS\n1 1,2\n2 3\n');
    });

    it('repeats an atom index when it carries weight > 1', () => {
        const a = { index: 4 };
        const collection = makeCollection([
            { name: 'B0', atoms: [a], expandedAtoms() { return [a, a]; } },
        ]);
        expect(generateBartender(collection)).toBe('BEADS\n1 5,5\n');
    });
});

describe('generatePyCGTOOL', () => {
    it('writes a comment header, section name, and one line per bead', () => {
        const collection = makeCollection(
            [{ name: 'R1', type: 'TC4', charge: 0, resname: 'NAPH',
               atoms: [{ index: 0 }, { index: 1 }],
               expandedAtoms: () => [{ index: 0 }, { index: 1 }] }],
            { 0: 'C8', 1: 'H8' },
        );
        const out = generatePyCGTOOL(collection);
        expect(out).toContain('[ NAPH ]');
        expect(out).toMatch(/R1\s+TC4\s+0\s+C8 H8/);
    });
});

describe('parseShakerMapping', () => {
    it('round-trips the real example mapping from example.ts', () => {
        const beads = parseShakerMapping(EXAMPLE_MAPPING);
        expect(beads).toHaveLength(6);

        const r11 = beads.find((b) => b.name === 'R11');
        expect(r11).toEqual({
            name: 'R11', type: 'SP2', charge: 0, atoms: ['C06', 'O08', 'N07', 'H13'],
        });

        const wr10 = beads.find((b) => b.name === 'WR10')!;
        expect(wr10.type).toBe('TN6a');
        expect(wr10.atoms).toEqual(['C0H', 'N0G', 'H0M']);
    });

    it('returns an empty array for text with no bead entries', () => {
        expect(parseShakerMapping('mapping = {}')).toEqual([]);
    });
});
