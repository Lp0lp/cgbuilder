import { describe, it, expect, beforeAll } from 'vitest';
import { Bead, BeadCollection } from '../scripts/bead.js';

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

function atom(index, x, y, z) {
    return { index, positionToVector3: () => ({ x, y, z }) };
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
});
