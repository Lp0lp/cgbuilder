import { describe, it, expect, beforeAll } from 'vitest';
import { Bead, BeadCollection } from '../scripts/bead.js';
import { buildStructure } from './helpers/mockStructure.js';

// Bead.center calls `new NGL.Vector3(...)`, referencing NGL as a browser
// global (loaded via <script> tag in index.html, not an ES import). Stand in
// with just enough of the real Vector3 API for the weighted-average math.
beforeAll(() => {
    globalThis.NGL = {
        Vector3: class {
            constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
            add(v) { this.x += v.x; this.y += v.y; this.z += v.z; return this; }
            divideScalar(s) { this.x /= s; this.y /= s; this.z /= s; return this; }
        },
    };
});

function atom(index, x, y, z, atomname = `A${index}`) {
    return { index, positionToVector3: () => ({ x, y, z }), atomname };
}

describe('Bead', () => {
    it('adding the same atom twice increases its weight without duplicating it', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        bead.addAtom(a);
        bead.addAtom(a);
        expect(bead.atoms).toHaveLength(1);
        expect(bead.atomWeights[0]).toBe(2);
    });

    it('shift-click (removeAtom) decrements weight, then removes at zero', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        bead.addAtom(a);
        bead.addAtom(a);
        bead.removeAtom(a);
        expect(bead.atoms).toHaveLength(1);
        expect(bead.atomWeights[0]).toBe(1);
        bead.removeAtom(a);
        expect(bead.atoms).toHaveLength(0);
        expect(bead.atomWeights[0]).toBeUndefined();
    });

    it('center is the unweighted average for equally-weighted atoms', () => {
        const bead = new Bead();
        bead.addAtom(atom(0, 0, 0, 0));
        bead.addAtom(atom(1, 10, 0, 0));
        const c = bead.center;
        expect(c.x).toBeCloseTo(5);
        expect(c.y).toBeCloseTo(0);
    });

    it('center pulls toward an atom weighted by repeated clicks', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        const b = atom(1, 10, 0, 0);
        bead.addAtom(a);
        bead.addAtom(b);
        bead.addAtom(b); // weight 2 on b -> center shifts toward b
        const c = bead.center;
        // weighted average: (0*1 + 10*2) / 3 = 6.667
        expect(c.x).toBeCloseTo(20 / 3);
    });

    it('charge coerces non-numeric input to 0', () => {
        const bead = new Bead();
        bead.charge = 'not-a-number';
        expect(bead.charge).toBe(0);
        bead.charge = '-1.5';
        expect(bead.charge).toBe(-1.5);
    });

    it('indexOf finds an atom by index, -1 if absent', () => {
        const bead = new Bead();
        const a = atom(3, 0, 0, 0);
        bead.addAtom(a);
        expect(bead.indexOf(a)).toBe(0);
        expect(bead.indexOf(atom(99, 0, 0, 0))).toBe(-1);
    });

    it('isAtomIn reflects membership', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        const b = atom(1, 0, 0, 0);
        bead.addAtom(a);
        expect(bead.isAtomIn(a)).toBe(true);
        expect(bead.isAtomIn(b)).toBe(false);
    });

    it('toggleAtom always adds/increments, never removes (see its own comment)', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        bead.toggleAtom(a);
        bead.toggleAtom(a);
        expect(bead.atoms).toHaveLength(1);
        expect(bead.atomWeights[0]).toBe(2);
    });

    it('expandedAtoms repeats each atom by its weight, default 1', () => {
        const bead = new Bead();
        const a = atom(0, 0, 0, 0);
        const b = atom(1, 0, 0, 0);
        bead.addAtom(a);
        bead.addAtom(b);
        bead.addAtom(b); // weight 2
        expect(bead.expandedAtoms()).toEqual([a, b, b]);
    });
});

describe('BeadCollection', () => {
    it('starts with one bead named B0', () => {
        const col = new BeadCollection();
        expect(col.beads).toHaveLength(1);
        expect(col.beads[0].name).toBe('B0');
        expect(col.currentBead).toBe(col.beads[0]);
    });

    it('names new beads with a strictly increasing index, even after removal', () => {
        const col = new BeadCollection();
        col.newBead(); // B1
        col.removeBead(1);
        const b2 = col.newBead();
        expect(b2.name).toBe('B2');
    });

    it('clearBeads resets state back to empty with no current bead', () => {
        const col = new BeadCollection();
        col.newBead();
        col.clearBeads();
        expect(col.beads).toHaveLength(0);
        expect(col.currentBead).toBeNull();
    });

    it('countBeadsForAtom counts how many beads reference a given atom', () => {
        const col = new BeadCollection();
        const a = atom(5, 0, 0, 0);
        col.beads[0].addAtom(a);
        const second = col.newBead();
        second.addAtom(a);
        expect(col.countBeadsForAtom(a)).toBe(2);
    });

    describe('atom naming', () => {
        it('atomName falls back to atom.atomname when no original names were set', () => {
            const col = new BeadCollection();
            expect(col.atomName(atom(0, 0, 0, 0, 'CA'))).toBe('CA');
        });

        it('atomName prefers the recorded original name, by index, once set', () => {
            const col = new BeadCollection();
            col.setOriginalAtomNames(['N1', 'C2', 'C3']);
            expect(col.atomName(atom(1, 0, 0, 0, 'fallback'))).toBe('C2');
            // Index 5 was never recorded -- falls back to atom.atomname.
            expect(col.atomName(atom(5, 0, 0, 0, 'C99'))).toBe('C99');
        });

        it('regression: recovers a PDB/GRO atom name NGL re-cased (e.g. "CL9" -> "Cl9")', () => {
            // setOriginalAtomNames/atomName exist specifically because NGL's
            // own parsed atom.atomname isn't guaranteed to match the source
            // file byte-for-byte -- it was re-casing a chlorine atom named
            // "CL9" in the source PDB to "Cl9", which broke anything that
            // needed the exact original name (e.g. round-tripping a Shaker
            // mapping, where the exported name must match what re-importing
            // expects). The fix reads names straight from the file's raw
            // text (see parsePDBAtomNames/parseGROAtomNames in
            // fileformats.js) and looks them up by index here instead of
            // trusting NGL's own casing.
            const col = new BeadCollection();
            col.setOriginalAtomNames(['C1', 'CL9']); // as read verbatim from the PDB
            const nglMisCasedAtom = atom(1, 0, 0, 0, 'Cl9'); // NGL's own (wrong) parse
            expect(col.atomName(nglMisCasedAtom)).toBe('CL9');
        });

        it('setOriginalAtomNames with a non-array input clears rather than throwing', () => {
            const col = new BeadCollection();
            col.setOriginalAtomNames(['N1']);
            col.setOriginalAtomNames(null);
            expect(col.atomName(atom(0, 0, 0, 0, 'fallback'))).toBe('fallback');
        });

        it('atomNames maps atomName over a list of atoms', () => {
            const col = new BeadCollection();
            col.setOriginalAtomNames(['N1', 'C2']);
            const names = col.atomNames([atom(0, 0, 0, 0), atom(1, 0, 0, 0)]);
            expect(names).toEqual(['N1', 'C2']);
        });

        it('expandedAtomNames repeats a weighted atom\'s name, using recorded original names', () => {
            const col = new BeadCollection();
            col.setOriginalAtomNames(['N1', 'C2']);
            const bead = col.beads[0];
            const a = atom(0, 0, 0, 0);
            const b = atom(1, 0, 0, 0);
            bead.addAtom(a);
            bead.addAtom(b);
            bead.addAtom(b); // weight 2
            expect(col.expandedAtomNames(bead)).toEqual(['N1', 'C2', 'C2']);
        });

        it('structureAtomNames returns atomName for every atom in structure order', () => {
            const col = new BeadCollection();
            col.setOriginalAtomNames(['N1', 'C2']);
            const structure = buildStructure([
                { element: 'N', atomname: 'fallback0' },
                { element: 'C', atomname: 'fallback1' },
            ]);
            expect(col.structureAtomNames(structure)).toEqual(['N1', 'C2']);
        });

        it('structureAtomNames returns an empty array for a structure with no eachAtom', () => {
            const col = new BeadCollection();
            expect(col.structureAtomNames(null)).toEqual([]);
            expect(col.structureAtomNames({})).toEqual([]);
        });
    });
});
