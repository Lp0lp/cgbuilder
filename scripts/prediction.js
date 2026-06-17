/* ===========================================================================
   Bead type prediction — AutoMartini M3 port
   ===========================================================================
   Fragment → delta_f (kJ/mol) lookup table — the full AutoMartini M3
   logP_smi_extended.dat (data/logP_smi_extended.dat), with duplicate keys
   resolved by FIRST occurrence (the original calibrated value; later duplicates
   are AutoMartini's dict-overwrite bug). Keys are kept verbatim in AutoMartini's
   notation: lowercase = aromatic atoms, uppercase = aliphatic; many are
   open-chain aromatic SMILES ("cc", "cn", "ncs") that RDKit cannot parse — these
   are matched directly against the aromatic-notation fragment key (see
   buildCanonTable / onPredictTypes). Do NOT "canonicalise" or hand-edit values. */
const FRAG_DELTA_F = {
    "CC":12.0, "CCC":14.2, "CCCC":18.9, "CC(C)(C)C":18.9,
    "C=C":6.4, "C=CC":8.4, "CC=C":8.4, "CCc":8.4,
    "cCC":8.4, "cc(c)C":13.4, "cc(c)-c":13.4, "CCCc":13.4,
    "CC(c)C":13.4, "ccCC":13.4, "C=CCC":13.4, "CC=CC":13.4,
    "CCC=C":13.4, "CC1(C)CC1":13.4, "cC1CC1":13.8, "cC(C)C":13.8,
    "C=C(C)C":13.8, "C=C(C)c":13.8, "CC(C)c":13.8, "CC(C)(C)c":13.8,
    "CC(C)C":14.8, "CCC(C)C":14.8, "cC":6.4, "Cc":6.4,
    "ccC":8.4, "cc":4.5, "Ccc":6.3,
    "ccc":6.3, "cc-c":6.3, "C=Cc":6.3, "C1CC1":9.2,
    "cCcc":11.2, "CC1CC1":11.2, "Cc(c)c":11.2, "cccc":11.2,
    "cccC":11.2, "C#C":3.6, "C#CC":5.3, "C#Cc":5.3,
    "CC#C":5.3, "C#CCC":10.1, "CC#CC":10.1, "CCC#C":10.1,
    "C#CC(C)C":10.1, "CS":3.6, "CCS":5.3, "CCCS":10.1,
    "cs":3.6, "cS":3.6, "ccs":5.3, "ccS":5.3,
    "cSc":5.3, "Scc":5.3, "Ccs":5.3, "c(CS)c":10.1,
    "C=C(C)S":10.1, "SC":3.6, "CSC":5.3, "SCC":5.3,
    "csc":5.3, "C=CS":5.3, "CSc":8.4, "CSCC":10.1,
    "cc(-c)s":10.1, "CCSC":10.1, "SCCC":10.1, "cc(C)s":10.1,
    "CC(C)S":13.4, "CC(C)(C)S":13.4, "CCF":2.7, "CCCF":4.3,
    "CCCCF":8.7, "FCC":2.7, "C(F)CC":4.3, "CC(F)C":4.3,
    "CC(F)CC":8.7, "cCl":5.4, "ccCl":8.0, "scCl":8.0,
    "CCl":5.4, "CCCl":8.0, "CCCCl":13.9, "ClCC":8.0,
    "C=CCl":8.0, "SCCl":8.0, "ocCl":4.3, "C=CF":9.4,
    "ncF":9.4, "C=C(c)Cl":13.9, "CC(Cl)C":13.9, "CC(C)Cl":13.9,
    "cc(c)Cl":13.9, "cc(s)Cl":13.9, "CC(F)F":13.9, "CC(C)(C)Cl":14.3,
    "CBr":5.2, "CCBr":7.2, "CCCBr":12.7, "cBr":5.2,
    "BrCC":7.2, "C=CBr":7.2, "ccBr":7.2, "ClCCl":7.2,
    "CC(Br)C":12.7, "C=C(C)Br":12.7, "C=C(C)Cl":12.7, "NC(F)F":12.7,
    "FC(F)S":12.7, "CC(C)(F)F":13.9, "cc(c)F":14.3, "CC(C)F":14.3,
    "OC(F)F":8.7, "CC(F)(F)F":8.7, "ClC1(Cl)CC1":8.7, "CC(Cl)Cl":8.7,
    "cC(F)(F)F":8.7, "OC(F)(F)F":8.7, "CC(O)F":8.7, "OCCl":4.3,
    "O=PCl":4.3, "FC(F)F":8.7, "I":7.6, "CF":7.6,
    "cF":7.6, "CI":7.6, "CCI":14.3, "IC":7.6,
    "cI":7.6, "ccI":9.4, "FCF":7.2, "ccF":8.0,
    "C(I)C":14.3, "Cn":2.3, "cN":-2.8, "CN":-4.1,
    "CCN":-2.5, "CC#N":0.2, "cC#N":1.1, "C=CN":-2.5,
    "CNC":3.6, "CCn":3.6, "CCCN":1.0, "CN(C)C":8.1,
    "CN(C)c":8.1, "CC(C)n":8.1, "CC(c)N":8.1, "CCC(C)(C)N":10.7,
    "NS":-4.1, "cn":-4.1, "C=N":-4.1,
    "O=s":-4.1, "NC=O":-6.9, "NCO":-5.4, "NCC":-0.9,
    "Ccn":6.0, "C#N":-1.2, "ccN":-2.5, "C=Nc":-2.5,
    "cNc":-2.5, "ccn":-1.0, "C=Cn":-1.0, "cnc":-1.0,
    "ncs":-1.0, "cns":-1.0, "NCS":-1.0, "cSC":5.3,
    "cCS":5.3, "C=NC":3.6, "cNC":3.6, "cno":1.1,
    "C=NO":1.1, "nco":1.1, "Nco":1.1, "cC=O":1.1,
    "con":1.1, "cOc":1.1, "nCO":-5.2, "Cnn":-5.2,
    "CNN":-5.2, "cNN":-5.2, "Ncs":-2.5, "O=CS":-2.5,
    "Cs=O":-2.5, "O=s=O":-3.7, "ncn":-3.7, "cnn":-3.7,
    "c-nn":-3.7, "nns":-3.7, "ncN":-6.9, "cnN":-6.9,
    "N=CN":-6.9, "C=Nn":-6.9, "CNc":-6.9, "c-cn":-2.5,
    "cnC":-2.5, "cn-c":-2.5, "CC=N":-2.5, "N=CS":-2.5,
    "ncS":-2.5, "csn":-2.5, "COn":-2.5, "C=NN":-1.0,
    "NCN":-5.4, "cCN":-4.2, "CNO":0.2, "nn":-6.1,
    "nN":-6.1, "Nn":-6.1, "NN":-7.8, "ns":-4.1,
    "no":-7.2, "N=O":-7.2, "O=P":-5.0, "NO":-5.0,
    "N=S=O":-1.0, "ccnn":0.2, "cncn":0.2, "nc(n)S":0.2,
    "Cc(n)n":0.2, "cc(n)n":0.2, "cn(-c)n":0.2, "nc(n)N":0.2,
    "CN(c)N":0.2, "C=C(N)N":0.2, "N[SH](=O)=O":-1.9, "CS(N)(=O)=O":-1.9,
    "N=C(N)N":-1.9, "c[N+](=O)[O-]":3.8, "CC(C)nn":3.8, "C[SH](=O)=O":0.2,
    "c[SH](=O)=O":0.2, "CS(C)(=O)=O":0.2, "N=[SH](N)=O":0.2, "NC(=O)S":0.2,
    "Cc(n)s":0.2, "c-c(n)s":0.2, "cc(N)s":0.2, "CS(C)=O":0.2,
    "cc(c)n":1.0, "CC(=N)S":1.0, "cc(n)S":1.0, "cc(-c)n":1.0,
    "C=C(C)n":1.0, "cn(c)C":1.0, "C#CCN":1.0, "CC(C)(C)n":8.1,
    "nc(N)s":-1.9, "cc(n)N":-1.9, "CC(=N)N":-1.9, "cn(C)n":-1.9,
    "N=C(N)S":-1.9, "c-c(n)n":-1.9, "ccC#N":4.3, "CC(C)C#N":4.3,
    "NC(N)=O":0.2, "C=C(N)O":0.2, "ncO":-5.2, "N=NN":-5.2,
    "CON":-1.0, "NC(=O)":-3.7, "N=CO":-3.7, "N[SH]=O":-3.7,
    "O=[N+][O-]":-1.8, "CC(C)=N":2.2, "cc(C)n":2.2, "CC(C)N":1.0,
    "CNC(C)C":1.0, "cN(C)C":1.0, "cC(C)N":1.0, "cccn":1.0,
    "ccnc":1.0, "nc(N)o":1.0, "CC(N)C":2.2, "C=C(C)N":2.2,
    "CC(N)=O":1.0, "CC(=N)O":1.0, "cC(C)=N":1.0, "cc(n)s":0.2,
    "cC(N)=O":0.2, "nc(n)O":0.2, "N=C(N)O":0.2, "C=C(N)n":0.2,
    "Cc(n)o":4.3, "cc(n)O":-1.9, "c-c(n)o":-1.9, "ccon":-1.1,
    "nC1CC1":8.1, "NC1CC1":8.1, "CCNC":8.1, "CCC(C)N":10.7,
    "cc(c)S":10.1, "cc(-c)o":7.8, "NCCO":7.8, "CC(C)(C)O":8.1,
    "CC(C)(C)N":8.1, "ccN(C)C":8.1, "CC(C)(c)O":8.1, "cc(c)N":3.8,
    "CCC#N":3.8, "CN(c)C":2.2, "ccNC":2.2, "CC(C)(N)O":4.3,
    "OC":-1.4, "CO":-1.4, "co":2.3, "COc":3.8,
    "COC":0.2, "CCO":0.2, "Cco":0.2, "cCO":-4.2,
    "cco":3.8, "CCOC":3.8, "COCC":3.8, "coc":3.8,
    "CP=O":-2.5, "C=O":-2.8, "CC=O":-1.0, "cC=N":-1.0,
    "CCC=O":2.2, "cC(O)O":2.2, "CC(C)O":9.9, "CC(c)O":-2.0,
    "cC(C)O":-2.0, "OCC":-1.0, "O=S=O":-3.7, "CC(=O)C":2.2,
    "C(=O)O":-1.2, "O=CO":-6.9, "CC(O)O":1.0, "CC(=O)O":-3.8,
    "cC(=O)O":-3.8, "ccOC":7.8, "COC(C)C":7.8, "cc(C)o":7.8,
    "C=C(C)O":7.8, "cc(c)o":7.8, "ccOc":7.8, "ccCO":4.3,
    "cOC":1.1, "CCCO":-2.0, "CC(O)C":-2.0, "cO":-6.1,
    "Ocs":-4.2, "ccO":-4.2, "C=CO":-4.2, "cc(c)O":-4.2,
    "cccO":-4.2,
};

/* Free energies of transfer (kJ/mol) for all Martini 3 bead types.
   Source: Souza et al. 2021 (https://doi.org/10.1038/s41592-021-01098-3). */
const DELTA_F = {
    C1:18.9, C2:14.8, C3:13.8, C4:13.4, C5:11.2, C6:10.1,
    N1:8.1,  N2:5.6,  N3:1.8,  N4:2.2,  N5:0.0,  N6:-1.1,
    P1:-2.0, P2:-3.8, P3:-5.1, P4:-7.4, P5:-9.1, P6:-9.2,
    X1:14.3, X2:12.7, X3:13.9, X4:8.7,
    N1d:10.7, N1a:10.7, N2d:7.8,  N2a:7.8,  N3d:3.8,  N3a:3.8,
    N4d:4.3,  N4a:4.3,  N5d:2.2,  N5a:2.2,  N6d:1.0,  N6a:1.0,
    P1d:0.2,  P1a:0.2,  P2d:-1.9, P2a:-1.9, P3d:-3.5, P3a:-3.5,
    P4d:-5.1, P4a:-5.1, P5d:-7.0, P5a:-7.0, P6d:-7.4, P6a:-7.4,
    Q1:-10.9, Q2:-15.1, Q3:-17.4, Q4:-18.8, Q5:-23.0, D:-26.8,

    SC1:14.2, SC2:9.9,  SC3:9.2,  SC4:8.4,  SC5:6.3,  SC6:5.3,
    SN1:3.6,  SN2:2.1,  SN3:-1.8, SN4:-0.9, SN5:-3.6, SN6:-4.2,
    SP1:-5.2, SP2:-6.9, SP3:-7.7, SP4:-9.8, SP5:-11.8, SP6:-12.0,
    SX1:9.4,  SX2:7.2,  SX3:8.0,  SX4:4.3,
    SN1d:6.0,  SN1a:6.0,  SN2d:3.8,  SN2a:3.8,  SN3d:0.2,  SN3a:0.2,
    SN4d:1.1,  SN4a:1.1,  SN5d:-1.0, SN5a:-1.0, SN6d:-2.5, SN6a:-2.5,
    SP1d:-3.7, SP1a:-3.7, SP2d:-5.4, SP2a:-5.4, SP3d:-6.1, SP3a:-6.1,
    SP4d:-7.8, SP4a:-7.8, SP5d:-9.5, SP5a:-9.5, SP6d:-9.6, SP6a:-9.6,
    SQ1:-10.6, SQ2:-14.3, SQ3:-18.0, SQ4:-18.2, SQ5:-18.2, SD:-36.4,

    TC1:12.0, TC2:7.8,  TC3:6.7,  TC4:6.4,  TC5:4.5,  TC6:3.6,
    TN1:2.3,  TN2:0.3,  TN3:-3.1, TN4:-2.9, TN5:-4.9, TN6:-6.1,
    TP1:-7.2, TP2:-8.8, TP3:-9.8, TP4:-12.1, TP5:-15.2, TP6:-14.8,
    TX1:7.6,  TX2:5.2,  TX3:5.4,  TX4:2.7,
    TN1d:3.9,  TN1a:3.9,  TN2d:2.3,  TN2a:2.3,  TN3d:-1.4, TN3a:-1.4,
    TN4d:-1.2, TN4a:-1.2, TN5d:-2.8, TN5a:-2.8, TN6d:-4.1, TN6a:-4.1,
    TP1d:-5.0, TP1a:-5.0, TP2d:-6.8, TP2a:-6.8, TP3d:-7.8, TP3a:-7.8,
    TP4d:-9.5, TP4a:-9.5, TP5d:-13.2, TP5a:-13.2, TP6d:-12.7, TP6a:-12.7,
    TQ1:-14.2, TQ2:-14.5, TQ3:-18.7, TQ4:-16.3, TQ5:-17.0, TD:-36.8,
};

// Build a canonicalized version of FRAG_DELTA_F using RDKit to normalize SMILES.
// Called once per predict session. Keys in the original table may be in aromatic
// notation (e.g. "cc") or non-canonical form; this ensures the lookup always
// succeeds when comparing against RDKit-canonicalized fragment SMILES.
export function buildCanonTable(RDKit) {
    const result = {};
    for (const [smi, df] of Object.entries(FRAG_DELTA_F)) {
        try {
            const mol = RDKit.get_mol(smi);
            if (!mol) {
                // Open-chain aromatic SMILES like "cc", "ccc" are invalid for
                // RDKit (no ring to close) but are legitimate AutoMartini table
                // keys. Keep them as-is for direct lookup by aromatic fragments.
                if (!(smi in result)) result[smi] = df;
                continue;
            }
            const canon = mol.get_smiles();
            mol.delete();
            if (canon && !(canon in result)) result[canon] = df;
        } catch (_) {
            // RDKit threw on an open-chain aromatic key — keep it verbatim so the
            // aromatic-notation fragment lookup can still match it.
            if (!(smi in result)) result[smi] = df;
        }
    }
    return result;
}

function _closestType(deltaF, candidates) {
    let best = candidates[0];
    let bestErr = Infinity;
    for (const t of candidates) {
        if (DELTA_F[t] === undefined) continue;
        const err = Math.abs(DELTA_F[t] - deltaF);
        if (err < bestErr) { bestErr = err; best = t; }
    }
    return best;
}

// Port of AutoMartini M3 determine_bead_type().
// deltaF: free energy of transfer (kJ/mol) — from table or Crippen conversion.
// heavyCount: non-H atoms in bead; inRing: any bead atom is aromatic.
export function determineBeadType({ deltaF, charge, hDonors, hAcceptors, hasHalogen, inRing, heavyCount }) {
    const sz = heavyCount <= 2 ? 'T' : heavyCount === 3 ? 'S' : '';
    const p = sz;

    let result;

    if (charge !== 0) {
        result = _closestType(deltaF,
            [`${p}Q1`,`${p}Q2`,`${p}Q3`,`${p}Q4`,`${p}Q5`,`${p}D`]);
    } else if (hAcceptors > 0 && hDonors === 0) {
        // Pure acceptor → 'a' label.
        result = _closestType(deltaF, p === 'T'
            ? ['TN1a','TN2a','TN3a','TN4a','TN5a','TN6a','TP1a','TP2a','TP3a','TP4a','TP5a','TP6a']
            : p === 'S'
            ? ['SN1a','SN2a','SN3a','SN4a','SN5a','SN6a','SP1a','SP2a','SP3a','SP4a','SP5a','SP6a']
            : ['N1a','N2a','N3a','N4a','N5a','N6a','P1a','P2a','P3a','P4a','P5a','P6a']);
    } else if (hDonors > 0 && hAcceptors === 0) {
        // Pure donor → 'd' label.
        result = _closestType(deltaF, p === 'T'
            ? ['TN1d','TN2d','TN3d','TN4d','TN5d','TN6d','TP1d','TP2d','TP3d','TP4d','TP5d','TP6d']
            : p === 'S'
            ? ['SN1d','SN2d','SN3d','SN4d','SN5d','SN6d','SP1d','SP2d','SP3d','SP4d','SP5d','SP6d']
            : ['N1d','N2d','N3d','N4d','N5d','N6d','P1d','P2d','P3d','P4d','P5d','P6d']);
    } else {
        // Either no H-bonding, OR both donor AND acceptor (a balanced group such
        // as an amide) → plain N/P/C series with no a/d label, matching Martini 3
        // where the a/d suffixes denote one-sided H-bonders only.
        let cands;
        if (p === 'T') {
            cands = ['TP6','TP5','TP4','TP3','TP2','TP1','TC6','TC5','TC4','TC3','TC2','TC1','TN6','TN5','TN4','TN3','TN2','TN1'];
            if (!inRing) cands = cands.filter(t => t !== 'TC5');
        } else if (p === 'S') {
            cands = ['SP6','SP5','SP4','SP3','SP2','SP1','SC6','SC5','SC4','SC3','SC2','SC1','SN6','SN5','SN4','SN3','SN2','SN1'];
        } else {
            cands = ['P6','P5','P4','P3','P2','P1','C6','C5','C4','C3','C2','C1','N6','N5','N4','N3','N2','N1'];
        }
        result = _closestType(deltaF, cands);
    }

    // Halogen overrides everything (same as AutoMartini)
    if (hasHalogen) {
        result = _closestType(deltaF, p === 'T' ? ['TX4','TX3','TX2','TX1']
                                    : p === 'S' ? ['SX4','SX3','SX2','SX1']
                                    :              ['X4','X3','X2','X1']);
    }

    return result;
}
