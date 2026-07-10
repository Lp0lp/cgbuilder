import { describe, it, expect } from 'vitest';
import { determineBeadType } from '../prediction.js';

describe('determineBeadType', () => {
    it('picks the exact C-series match when deltaF lines up with a known value', () => {
        // DELTA_F.C1 = 18.9 exactly -> zero error, must win regardless of order.
        const result = determineBeadType({
            deltaF: 18.9, charge: 0, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).toBe('C1');
    });

    it('routes pure acceptors (HBA>0, HBD=0) to the "a"-suffixed series', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 0, hAcceptors: 1,
            hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).toMatch(/a$/);
    });

    it('routes pure donors (HBD>0, HBA=0) to the "d"-suffixed series', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 1, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).toMatch(/d$/);
    });

    it('does not suffix a/d when both donor and acceptor are present (e.g. amide)', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 0, hDonors: 1, hAcceptors: 1,
            hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).not.toMatch(/[ad]$/);
    });

    it('routes any non-zero charge to the Q/D series', () => {
        for (const charge of [1, -1, 2, -2]) {
            const result = determineBeadType({
                deltaF: -5.0, charge, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
            });
            expect(result).toMatch(/^(Q[1-5]|D)$/);
        }
    });

    it('breaks the SQ4/SQ5 tie (both -18.2 in DELTA_F) in favour of SQ5', () => {
        // Regression: a carboxylate-like ion-table entry landing on S-size at
        // exactly -18.2 used to resolve to SQ4 instead of the intended SQ5,
        // since they're the only exact tie anywhere in DELTA_F.
        const result = determineBeadType({
            deltaF: -18.2, charge: -1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 3, ringOrBranched: false,
        });
        expect(result).toBe('SQ5');
    });

    describe('divalent+ charge defaults to D', () => {
        it('uses the Q1-Q5 ladder for monovalent charge (+1/-1)', () => {
            for (const charge of [1, -1]) {
                const result = determineBeadType({
                    deltaF: -5.0, charge, hDonors: 0, hAcceptors: 0,
                    hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
                });
                expect(result).toMatch(/^Q[1-5]$/);
            }
        });

        it('bypasses deltaF entirely and returns D for |charge| >= 2', () => {
            for (const charge of [2, -2, 3, -3]) {
                // deltaF deliberately set to something that would clearly win
                // a Q-tier match if the ladder were searched, to prove the D
                // override short-circuits before any matching happens.
                const result = determineBeadType({
                    deltaF: -10.9, charge, hDonors: 0, hAcceptors: 0,
                    hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
                });
                expect(result).toBe('D');
            }
        });

        it('applies the size prefix to D the same way as any other bead', () => {
            const tiny = determineBeadType({
                deltaF: -5.0, charge: 2, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 1, ringOrBranched: false,
            });
            expect(tiny).toBe('TD');

            const small = determineBeadType({
                deltaF: -5.0, charge: -2, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 3, ringOrBranched: false,
            });
            expect(small).toBe('SD');
        });
    });

    it('halogen presence overrides everything else, including charge', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 1, hAcceptors: 1,
            hasHalogen: true, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).toMatch(/^X[1-4]$/);
    });

    it('uses the T- (tiny) prefix for a weighted count of 2', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 2, ringOrBranched: false,
        });
        expect(result).toMatch(/^T(Q[1-5]|D)$/);
    });

    it('uses the S- (small) prefix for a weighted count of 3', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 3, ringOrBranched: false,
        });
        expect(result).toMatch(/^S(Q[1-5]|D)$/);
    });

    it('uses no size prefix (regular) for a weighted count of 4 when linear (not ring/branched)', () => {
        const result = determineBeadType({
            deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
        });
        expect(result).toMatch(/^(Q[1-5]|D)$/);
    });

    it('excludes TC5 (non-ring-only tiny carbon) when the fragment is not in a ring', () => {
        // TC5 is filtered out of the tiny candidate list unless inRing=true
        // (see prediction.ts comment on the non-bonding branch).
        const result = determineBeadType({
            deltaF: 4.5, charge: 0, hDonors: 0, hAcceptors: 0,
            hasHalogen: false, inRing: false, weightedHeavyCount: 2, ringOrBranched: false,
        });
        expect(result).not.toBe('TC5');
    });

    describe('ring/branched bead-size override', () => {
        it('downgrades a weighted count of 4 from R to S when ring or branched', () => {
            const result = determineBeadType({
                deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: true,
            });
            expect(result).toMatch(/^S(Q[1-5]|D)$/);
        });

        it('does not apply the override for a weighted count of 4 when linear', () => {
            const result = determineBeadType({
                deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 4, ringOrBranched: false,
            });
            expect(result).toMatch(/^(Q[1-5]|D)$/); // no S/T prefix
        });

        it('leaves a weighted count of 5 as R even when ring or branched (override only applies at exactly 4)', () => {
            const result = determineBeadType({
                deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 5, ringOrBranched: true,
            });
            expect(result).toMatch(/^(Q[1-5]|D)$/); // no S/T prefix
        });

        it('leaves weighted counts of 2 and 3 unaffected by ring/branched (already T/S regardless)', () => {
            const tiny = determineBeadType({
                deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 2, ringOrBranched: true,
            });
            expect(tiny).toMatch(/^T(Q[1-5]|D)$/);

            const small = determineBeadType({
                deltaF: -5.0, charge: 1, hDonors: 0, hAcceptors: 0,
                hasHalogen: false, inRing: false, weightedHeavyCount: 3, ringOrBranched: true,
            });
            expect(small).toMatch(/^S(Q[1-5]|D)$/);
        });
    });
});
