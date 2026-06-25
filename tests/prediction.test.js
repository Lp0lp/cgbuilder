import { describe, it, expect } from 'vitest';
import { determineBeadType } from '../scripts/prediction.js';

describe('determineBeadType', () => {
    it('picks the exact C-series match when deltaF lines up with a known value', () => {
        // DELTA_F.C1 = 18.9 exactly -> zero error, must win regardless of order.
        const result = determineBeadType({
            deltaF: 18.9, charge: 0, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 4,
        });
        expect(result).toBe('C1');
    });

    it('routes pure acceptors (HBA>0, HBD=0) to the "a"-suffixed series', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 0, hAcceptors: 1,
            hasHalogen: false, inRing: false, heavyCount: 4,
        });
        expect(result).toMatch(/a$/);
    });

    it('routes pure donors (HBD>0, HBA=0) to the "d"-suffixed series', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 1, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 4,
        });
        expect(result).toMatch(/d$/);
    });

    it('does not suffix a/d when both donor and acceptor are present (e.g. amide)', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 1, hAcceptors: 1,
            hasHalogen: false, inRing: false, heavyCount: 4,
        });
        expect(result).not.toMatch(/[ad]$/);
    });

    it('routes any non-zero charge to the Q/D series', () => {
        for (const charge of [1, -1, 2, -2]) {
            const result = determineBeadType({
                deltaF: -5.0, charge, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, heavyCount: 4,
            });
            expect(result).toMatch(/^(Q[1-5]|D)$/);
        }
    });

    it('halogen presence overrides everything else, including charge', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 1, hAcceptors: 1,
            hasHalogen: true, inRing: false, heavyCount: 4,
        });
        expect(result).toMatch(/^X[1-4]$/);
    });

    it('uses the T- (tiny) prefix for 2-heavy-atom fragments', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 2,
        });
        expect(result).toMatch(/^T(Q[1-5]|D)$/);
    });

    it('uses the S- (small) prefix for 3-heavy-atom fragments', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 3,
        });
        expect(result).toMatch(/^S(Q[1-5]|D)$/);
    });

    it('uses no size prefix (regular) for 4+ heavy-atom fragments', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 4,
        });
        expect(result).toMatch(/^(Q[1-5]|D)$/);
    });

    it('excludes TC5 (non-ring-only tiny carbon) when the fragment is not in a ring', () => {
        // TC5 is filtered out of the tiny candidate list unless inRing=true
        // (see prediction.js comment on the non-bonding branch).
        const result = determineBeadType({
            deltaF: 4.5, charge: 0, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, heavyCount: 2,
        });
        expect(result).not.toBe('TC5');
    });
});
