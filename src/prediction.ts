/* ===========================================================================
   Bead type prediction — inspired by AutoMartini M3
   ===========================================================================
   Maps one fragment's physicochemical properties to a specific Martini 3
   bead type code (e.g. "SP2a"), largely inspired by AutoMartini's own
   determine_bead_type() logic (Szczuka et al. 2025,
   https://doi.org/10.1021/acs.jctc.5c01178; reference implementation
   https://github.com/Martini-Force-Field-Initiative/Automartini_M3/tree/main ):

     1. A free energy of transfer (deltaF, kJ/mol) for the fragment is
        looked up in FRAG_DELTA_F/ION_DELTA_F (via buildCanonTable's
        RDKit-canonicalized lookup), or falls back to a Crippen-logP-derived
        estimate when no table entry matches (computed by the caller, not
        in this file).
     2. determineBeadType picks a size prefix (T/S/none="regular") from the
        fragment's weighted heavy-atom count, with a ring/branch downgrade
        rule (see its own comment).
     3. It picks a candidate list of bead-type codes from charge and H-bond
        donor/acceptor pattern: divalent+ charge -> D, monovalent charge ->
        the Q1-Q5 ladder, pure acceptor/donor -> the "a"/"d"-suffixed
        series, anything else -> the plain C/N/P series.
     4. _closestType searches that candidate list against DELTA_F (the
        literature-calibrated reference free energy for every real bead
        type) and returns whichever candidate's own deltaF is numerically
        closest to the fragment's.
     5. A halogen-containing fragment overrides the result with the X1-X4
        series regardless of any of the above, matching AutoMartini.

   Fragment → delta_f (kJ/mol) lookup table — the full AutoMartini M3
   logP_smi_extended.dat. Keys are kept verbatim in AutoMartini's
   notation: lowercase = aromatic atoms, uppercase = aliphatic; many are
   open-chain aromatic SMILES ("cc", "cn", "ncs") that RDKit cannot parse — these
   are matched directly against the aromatic-notation fragment key (see
   buildCanonTable / onPredictTypes). */
import type { BeadTypeProps, RDKitModule } from './types.js';

const FRAG_DELTA_F: Record<string, number> = {
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

/* Common charged ions/functional groups — NOT from AutoMartini's table, our
   own additions, kept separate from FRAG_DELTA_F above so the verbatim
   AutoMartini data stays traceable.

   deltaF computed via Crippen logP collapses charged fragments to the same
   Q1/SQ1/TQ1 tier almost regardless of group (it's not calibrated for ionic
   species at all — verified empirically). These entries sidestep that:
   each deltaF value is chosen to land EXACTLY on the reference value of the
   plain Q-tier (size prefix handled separately, same dynamic logic as every
   other bead) that the Martini 3 SI lists for that ion/group, per its
   "Ions" table. Several of these are genuinely n/p-labelled in the SI
   (e.g. carboxylate -> Q5n, ammonium -> Q5p); the polarity label is
   dropped here since DELTA_F has no n/p reference values and
   determineBeadType doesn't search n/p candidates — landing on the right
   plain tier is a closer prediction than collapsing to Q1, and the user can
   add the label by hand. Chain-length variants are enumerated explicitly
   (not pattern-matched) since Martini beads only ever cover ~2-5 heavy
   atoms, a small, fully enumerable space. */
const ION_DELTA_F: Record<string, number> = {
    "[NH4+]":-17.0,
    "C[NH3+]":-16.3, "CC[NH3+]":-18.2, "CCC[NH3+]":-18.8, "CCCC[NH3+]":-18.8,
    "C[NH2+]C":-18.0, "CC[NH2+]C":-17.4, "CCC[NH2+]C":-17.4,
    "C[NH+](C)C":-14.3, "CC[NH+](C)C":-15.1,
    "C[N+](C)(C)C":-15.1, "C[N+](C)(C)CC":-10.9,
    "C[P+](C)(C)C":-10.9,
    "NC(=[NH2+])N":-18.0,
    "[O-]C=O":-18.2, "CC(=O)[O-]":-18.2, "CCC(=O)[O-]":-23.0,
    "CS(=O)(=O)[O-]":-18.8,
    "[O-]P(=O)(O)O":-23.0,
    "[B-](F)(F)(F)F":-15.1,
    "[P-](F)(F)(F)(F)(F)F":-10.9,
    "[S-]C#N":-10.6,
    "[N+](=O)([O-])[O-]":-18.0,
    "[O-]Cl(=O)(=O)=O":-15.1,
};

/* Reference free energy of transfer (kJ/mol) for every standard Martini 3
   bead type, keyed by its type code (e.g. "SP2a"). This is the table
   _closestType searches to turn a fragment's own deltaF into the
   nearest-matching real bead type. Source: Souza et al. 2021
   (https://doi.org/10.1038/s41592-021-01098-3).
   Note: SQ4 and SQ5 share the exact same value (-18.2) — the only exact
   tie anywhere in this table; see determineBeadType's Q-series ordering
   comment for how that tie is resolved. */
const DELTA_F: Record<string, number> = {
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

/**
 * Canonicalize every key in FRAG_DELTA_F/ION_DELTA_F through RDKit, so a
 * fragment's own RDKit-canonical SMILES can be looked up directly instead
 * of needing to match the table's original (often non-canonical) notation.
 * Open-chain aromatic keys ("cc", "cn", ...) that RDKit can't parse at all
 * are kept verbatim, for direct string-matching against the aromatic-
 * notation fragment key instead (see fragmentToSmiles's aromaticNotation
 * option in chemistry.ts). Called once per predict session, not per bead.
 * @param RDKit - the loaded RDKit_minimal module (see rdkit.ts)
 * @returns canonical (or verbatim, for aromatic-notation keys) SMILES ->
 *   deltaF (kJ/mol)
 */
export function buildCanonTable(RDKit: RDKitModule): Record<string, number> {
    const result: Record<string, number> = {};
    for (const table of [FRAG_DELTA_F, ION_DELTA_F]) {
        for (const [smi, df] of Object.entries(table)) {
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
    }
    return result;
}

/**
 * The candidate bead type (from `candidates`) whose own reference deltaF
 * (looked up in DELTA_F) is numerically closest to the fragment's deltaF.
 * On an exact tie, the FIRST candidate encountered wins — determineBeadType
 * relies on this for its Q5-before-Q4 ordering (see its own comment).
 * @param deltaF - the fragment's free energy of transfer (kJ/mol)
 * @param candidates - bead type codes to search, in order
 * @returns the closest-matching candidate
 */
function _closestType(deltaF: number, candidates: string[]): string {
    let best = candidates[0];
    let bestErr = Infinity;
    for (const t of candidates) {
        if (DELTA_F[t] === undefined) continue;
        const err = Math.abs(DELTA_F[t] - deltaF);
        if (err < bestErr) { bestErr = err; best = t; }
    }
    return best;
}

/**
 * Apply a Martini size prefix (T/S/'') and H-bond suffix (a/d/'') to every
 * bare code in `codes`, e.g. _series('S', ['N1','N2'], 'a') -> ['SN1a','SN2a'].
 * Lets determineBeadType's per-size-class candidate lists be written once
 * instead of by hand for every T/S/plain combination.
 * @param prefix - size class, 'T'/'S'/''
 * @param codes - bare bead-type codes, e.g. 'N1'
 * @param suffix - H-bond label, 'a'/'d'/''
 */
function _series(prefix: string, codes: string[], suffix = ''): string[] {
    return codes.map((c) => `${prefix}${c}${suffix}`);
}

const _NP_CODES = ['N1','N2','N3','N4','N5','N6','P1','P2','P3','P4','P5','P6'];
const _PCN_CODES = ['P6','P5','P4','P3','P2','P1','C6','C5','C4','C3','C2','C1','N6','N5','N4','N3','N2','N1'];
const _X_CODES = ['X4','X3','X2','X1'];

/**
 * Port of AutoMartini M3's determine_bead_type(): map one fragment's
 * physicochemical properties to a specific Martini 3 bead type code. Picks
 * a size prefix from weightedHeavyCount/ringOrBranched, then a candidate
 * list of bead-type codes from charge and H-bond donor/acceptor pattern
 * (see the module-level comment for the full decision order), and returns
 * whichever candidate's own reference deltaF (in DELTA_F) is closest to
 * this fragment's.
 * @param props - see BeadTypeProps
 * @returns the predicted bead type code, e.g. "SP2a"
 */
export function determineBeadType({
    deltaF, charge, hDonors, hAcceptors, hasHalogen, inRing, weightedHeavyCount, ringOrBranched,
}: BeadTypeProps): string {
    let sz = weightedHeavyCount <= 2 ? 'T' : weightedHeavyCount === 3 ? 'S' : '';
    if (sz === '' && weightedHeavyCount === 4 && ringOrBranched) sz = 'S';
    const p = sz;

    let result: string;

    if (Math.abs(charge) >= 2) {
        // Divalent+ ions (Mg2+, Ca2+, phosphate2-, ...) are D-beads by Martini
        // 3 convention regardless of deltaF -- the Q1-Q5 ladder is for
        // monovalent charge only, and deltaF-based matching is unreliable
        // for charged species anyway (verified empirically elsewhere).
        result = `${p}D`;
    } else if (charge !== 0) {
        // Q5 before Q4: SQ4 and SQ5 are the only exact tie anywhere in
        // DELTA_F (both -18.2). _closestType keeps whichever candidate it
        // sees first on a tie, so this ordering picks Q5 there instead of
        // Q4 -- harmless everywhere else, since no other deltaF value ties.
        result = _closestType(deltaF, _series(p, ['Q1', 'Q2', 'Q3', 'Q5', 'Q4', 'D']));
    } else if (hAcceptors > 0 && hDonors === 0) {
        // Pure acceptor → 'a' label.
        result = _closestType(deltaF, _series(p, _NP_CODES, 'a'));
    } else if (hDonors > 0 && hAcceptors === 0) {
        // Pure donor → 'd' label.
        result = _closestType(deltaF, _series(p, _NP_CODES, 'd'));
    } else {
        // Either no H-bonding, OR both donor AND acceptor (a balanced group such
        // as an amide) → plain N/P/C series with no a/d label, matching Martini 3
        // where the a/d suffixes denote one-sided H-bonders only.
        let cands = _series(p, _PCN_CODES);
        if (p === 'T' && !inRing) cands = cands.filter((t) => t !== 'TC5');
        result = _closestType(deltaF, cands);
    }

    // Halogen overrides everything (same as AutoMartini)
    if (hasHalogen) {
        result = _closestType(deltaF, _series(p, _X_CODES));
    }

    return result;
}
