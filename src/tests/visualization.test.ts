import { describe, it, expect, vi } from 'vitest';

// visualization.ts runtime-imports the NGL adapter (../ngl.js -> the `ngl`
// package) and rdkit.ts (which imports the RDKit WASM asset) at module load.
// Neither is needed to test the two pure functions below, and both are
// browser/bundler-only, so stub them out. vitest hoists these above the
// import below, so the real modules are never evaluated.
vi.mock('../ngl.js', () => ({ NGL: {} }));
vi.mock('../rdkit.js', () => ({ loadRDKit: () => Promise.reject(new Error('RDKit not needed in this test')) }));

import { typeColor, findParentWithClass } from '../visualization.js';

// Minimal duck-typed stand-in for a DOM Element, just enough surface for
// findParentWithClass's classList.contains/parentElement walk -- this
// project's tests run in vitest's node environment with no DOM, so a real
// Element isn't available (consistent with mockStructure.ts's approach for
// NGL elsewhere in this test suite).
function fakeElement(classes: string[], parent: Element | null = null): Element {
    const set = new Set(classes);
    return { classList: { contains: (c: string) => set.has(c) }, parentElement: parent } as unknown as Element;
}

describe('typeColor', () => {
    it('reads the polarity class from a plain (no size-prefix) type code', () => {
        expect(typeColor('C1')).toEqual([0.55, 0.55, 0.55]);
        expect(typeColor('N3')).toEqual([0.29, 0.56, 0.85]);
        expect(typeColor('P4')).toEqual([0.91, 0.30, 0.24]);
        expect(typeColor('Q2')).toEqual([0.95, 0.61, 0.07]);
        expect(typeColor('D')).toEqual([0.95, 0.61, 0.07]);
        expect(typeColor('X1')).toEqual([0.18, 0.80, 0.44]);
    });

    it('reads the polarity class from the SECOND letter for size-prefixed codes', () => {
        expect(typeColor('SP2a')).toEqual([0.91, 0.30, 0.24]); // S + P
        expect(typeColor('TC5')).toEqual([0.55, 0.55, 0.55]);  // T + C
        expect(typeColor('TQ1')).toEqual([0.95, 0.61, 0.07]);  // T + Q
    });

    it('is case-insensitive', () => {
        expect(typeColor('sp2a')).toEqual(typeColor('SP2a'));
    });

    it('falls back to light grey for an unknown or placeholder type', () => {
        expect(typeColor('TYPe')).toEqual([0.75, 0.75, 0.75]);
        expect(typeColor('')).toEqual([0.75, 0.75, 0.75]);
        expect(typeColor(undefined as unknown as string)).toEqual([0.75, 0.75, 0.75]);
    });
});

describe('findParentWithClass', () => {
    it('returns the element itself when it already has the class', () => {
        const el = fakeElement(['bead-view']);
        expect(findParentWithClass(el, 'bead-view')).toBe(el);
    });

    it('walks up through ancestors to find the nearest match', () => {
        const card = fakeElement(['bead-view']);
        const wrap = fakeElement(['field'], card);
        const input = fakeElement(['bead-name'], wrap);
        expect(findParentWithClass(input, 'bead-view')).toBe(card);
    });

    it('returns null when no ancestor (including the root) has the class', () => {
        const root = fakeElement(['page'], null);
        const child = fakeElement(['field'], root);
        expect(findParentWithClass(child, 'bead-view')).toBeNull();
    });
});
