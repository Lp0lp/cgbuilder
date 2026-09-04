import { describe, it, expect } from 'vitest';
import {
  shrakeRupley,
  beadSizeClass,
  beadRadius,
  BEAD_RADII,
} from '../sasa.js';
import type { Bead } from '../bead.js';

// beadRadius only reads `.type`; stand in a minimal object for a real Bead.
const bead = (type: string) => ({ type }) as unknown as Bead;

describe('beadSizeClass', () => {
  it('reads the size letter from the bead type prefix', () => {
    expect(beadSizeClass('SP2')).toBe('S');
    expect(beadSizeClass('TC5')).toBe('T');
    expect(beadSizeClass('UNK')).toBe('U');
  });

  it('defaults to regular (R) for types with no S/T/U prefix', () => {
    expect(beadSizeClass('P2')).toBe('R');
    expect(beadSizeClass('Q1')).toBe('R');
  });

  it('defaults to regular (R) for empty/missing types', () => {
    expect(beadSizeClass(null as unknown as string)).toBe('R');
    expect(beadSizeClass('')).toBe('R');
  });
});

describe('beadRadius', () => {
  it('maps each size class to its Martini vdW radius', () => {
    expect(beadRadius(bead('TC5'))).toBe(BEAD_RADII.T);
    expect(beadRadius(bead('SP2'))).toBe(BEAD_RADII.S);
    expect(beadRadius(bead('P2'))).toBe(BEAD_RADII.R);
    expect(beadRadius(bead('UNK'))).toBe(0);
  });
});

describe('shrakeRupley', () => {
  it('gives a fully isolated sphere its full theoretical surface area', () => {
    const r = 1.75,
      probe = 1.91;
    const shellR = r + probe;
    const expected = 4 * Math.PI * shellR * shellR;
    const result = shrakeRupley([[0, 0, 0, r]], probe, 2000);
    expect(result).toBeCloseTo(expected, 0);
  });

  it('two far-apart spheres sum to roughly twice the isolated-sphere area', () => {
    const r = 1.75,
      probe = 1.91;
    const isolated = shrakeRupley([[0, 0, 0, r]], probe, 2000);
    const farApart = shrakeRupley(
      [
        [0, 0, 0, r],
        [1000, 0, 0, r],
      ],
      probe,
      2000,
    );
    expect(farApart).toBeCloseTo(2 * isolated, -1);
  });

  it('overlapping spheres expose strictly less area than far-apart spheres', () => {
    const r = 1.75,
      probe = 1.91;
    const overlapping = shrakeRupley(
      [
        [0, 0, 0, r],
        [1.0, 0, 0, r],
      ],
      probe,
      2000,
    );
    const farApart = shrakeRupley(
      [
        [0, 0, 0, r],
        [1000, 0, 0, r],
      ],
      probe,
      2000,
    );
    expect(overlapping).toBeLessThan(farApart);
  });

  it('returns 0 for an empty particle list', () => {
    expect(shrakeRupley([], 1.91)).toBe(0);
  });
});
