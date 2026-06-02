const {
    Bead,
    BeadCollection,
    generateNDX,
    generateMap,
    generatePythonAssignments,
    generateGRO,
    parsePDBAtomNames,
    parseGROAtomNames,
    parseOriginalAtomNames,
    bondAwareRepresentationParams,
} = require('../scripts/main.js');

describe('Bead', () => {
    it('tracks atom weights and removes atoms when weight reaches zero', () => {
        const bead = new Bead();
        const atom = {
            index: 7,
            resname: 'MOL',
            resno: 3,
            atomname: 'C1',
            positionToVector3() {
                return { add() {}, divideScalar() {} };
            },
        };

        bead.addAtom(atom);
        bead.addAtom(atom);
        expect(bead.atoms).toHaveLength(1);
        expect(bead.atomWeights[7]).toBe(2);

        bead.removeAtom(atom);
        expect(bead.atoms).toHaveLength(1);
        expect(bead.atomWeights[7]).toBe(1);

        bead.removeAtom(atom);
        expect(bead.atoms).toHaveLength(0);
        expect(bead.atomWeights[7]).toBeUndefined();
    });
});

describe('BeadCollection atom naming', () => {
    it('uses original atom names when available and falls back otherwise', () => {
        const collection = new BeadCollection();
        const a0 = { index: 0, atomname: 'FALLBACK0' };
        const a1 = { index: 1, atomname: 'FALLBACK1' };

        collection.setOriginalAtomNames(['ORIG0']);

        expect(collection.atomName(a0)).toBe('ORIG0');
        expect(collection.atomName(a1)).toBe('FALLBACK1');
    });

    it('expands atom names according to bead weights', () => {
        const collection = new BeadCollection();
        const atom = { index: 0, atomname: 'C1' };
        const bead = collection.currentBead;

        bead.addAtom(atom);
        bead.addAtom(atom);
        collection.setOriginalAtomNames(['CA']);

        expect(collection.expandedAtomNames(bead)).toEqual(['CA', 'CA']);
    });
});

describe('Export generators', () => {
    it('builds NDX output', () => {
        const bead = { name: 'B0', atoms: [{ index: 0 }, { index: 3 }] };
        const out = generateNDX({ beads: [bead] });
        expect(out).toContain('[ B0 ]');
        expect(out).toContain('1 4');
    });

    it('builds MAP output using atom names from collection', () => {
        const atomA = { index: 1, atomname: 'A' };
        const atomB = { index: 0, atomname: 'B' };
        const collection = {
            beads: [{ name: 'B0', atoms: [atomA, atomB] }],
            atomName(atom) {
                return atom.atomname;
            },
        };

        const out = generateMap(collection);
        expect(out).toContain('[ atoms ]');
        expect(out).toContain('A');
        expect(out).toContain('B');
        expect(out).toContain('B0');
    });

    it('builds Python assignment output with duplicated mapped atoms', () => {
        const bead = {
            name: 'B0',
            type: 'C1',
            charge: -0.2,
            atoms: [{ index: 0, resname: 'MOL' }],
            atomWeights: { 0: 2 },
            resname: 'MOL',
        };
        const collection = {
            beads: [bead],
            expandedAtomNames() {
                return ['CA', 'CA'];
            },
        };

        const out = generatePythonAssignments(collection);
        expect(out).toContain('"MOL"');
        expect(out).toContain('"B0"');
        expect(out).toContain("'CA','CA'");
    });

    it('builds GRO output from bead center values', () => {
        const bead = {
            resid: 1,
            resname: 'MOL',
            name: 'B0',
            center: { x: 10, y: 20, z: 30 },
        };

        const out = generateGRO({ beads: [bead] });
        expect(out).toContain('Generated with cgbuilder');
        expect(out).toContain('    1');
        expect(out).toContain('   1.000');
        expect(out).toContain('   2.000');
        expect(out).toContain('   3.000');
    });
});

describe('Parsing helpers', () => {
    it('parses atom names from PDB content', () => {
        const pdb = [
            'ATOM      1  CA  ALA A   1      11.104  13.207  14.099  1.00 20.00           C',
            'HETATM    2  O1  LIG A   1      12.000  10.000  11.000  1.00 20.00           O',
        ].join('\n');

        expect(parsePDBAtomNames(pdb)).toEqual(['CA', 'O1']);
    });

    it('parses atom names from GRO content', () => {
        const gro = [
            'Test',
            '2',
            '    1MOL     CA    1   0.000   0.000   0.000',
            '    1MOL     CB    2   0.000   0.000   0.000',
            '   1.00000   1.00000   1.00000',
        ].join('\n');

        expect(parseGROAtomNames(gro)).toEqual(['CA', 'CB']);
    });

    it('dispatches parser based on file extension', () => {
        const pdb =
            'ATOM      1  CA  ALA A   1      11.104  13.207  14.099  1.00 20.00           C\n';
        const gro = [
            'Test',
            '1',
            '    1MOL     CA    1   0.000   0.000   0.000',
            '   1.00000   1.00000   1.00000',
        ].join('\n');

        expect(parseOriginalAtomNames(pdb, 'foo.pdb')).toEqual(['CA']);
        expect(parseOriginalAtomNames(gro, 'foo.gro')).toEqual(['CA']);
        expect(parseOriginalAtomNames('x', 'foo.sdf')).toEqual([]);
    });
});

describe('Representation params helper', () => {
    it('returns defaults and applies overrides', () => {
        const out = bondAwareRepresentationParams({
            bondScale: 0.9,
            sele: 'all',
        });
        expect(out.multipleBond).toBe(true);
        expect(out.bondSpacing).toBe(1);
        expect(out.bondScale).toBe(0.9);
        expect(out.sele).toBe('all');
    });
});
