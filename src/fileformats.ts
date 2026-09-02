/* ===========================================================================
   File I/O — export generators, import parsers, small browser helpers
   ===========================================================================
   Three unrelated concerns share this file because each is small enough not
   to need its own module:
     1. Output generators (generateNDX, generateMap, generatePythonAssignments,
        generateGRO, generateBartender) — serialise the current bead mapping
        into the formats downstream tools expect (GROMACS .ndx/.gro, a Martini
        .map, a Shaker Python assignments dict, a Bartender mapping file).
     2. Import parsers (parsePDBAtomNames/parseGROAtomNames/
        parseOriginalAtomNames/readOriginalAtomNames, parseShakerMapping) —
        the reverse direction: recovering the original AA atom names from a
        just-loaded file (NGL's own parsed atom names aren't guaranteed to
        match the source file byte-for-byte), and re-importing a previously
        exported Shaker mapping.
     3. Small browser utility helpers (download, copyTextToClipboard) and one
        NGL representation-params helper (bondAwareRepresentationParams) that
        don't have a more specific home elsewhere. */
import type { BeadCollection } from './bead.js';
import type { AtomProxy, BeadDef, RepresentationParams } from './types.js';

/**
 * GROMACS .ndx text: one `[ bead name ]` group per bead, listing its atoms'
 * 1-based (GROMACS convention) indices.
 * @param collection - BeadCollection
 */
export function generateNDX(collection: BeadCollection): string {
    let ndx = "";
    for (const bead of collection.beads) {
        ndx += "[ " + bead.name + " ]\n";
        for (const atom of bead.atoms) {
            ndx += (atom.index + 1) + " ";
        }
        ndx += "\n\n";
    }
    return ndx;
}

/**
 * Martini .map text: a `[ to ]`/`[ martini ]` header followed by an
 * `[ atoms ]` section, one line per distinct AA atom NAME listing which
 * bead(s) it contributes to. Grouping is by atom name rather than atom
 * index — this assumes atom names are unique within the structure (the
 * same assumption the app's own atom-name-uniqueness warning exists to
 * flag; if two different atoms share a name, this would wrongly merge
 * them into one .map line). Atom order in the output follows the original
 * structure's atom index, not bead order.
 * @param collection - BeadCollection
 */
export function generateMap(collection: BeadCollection): string {
    let output = "[ to ]\nmartini\n\n[ martini ]\n";
    const atomToBeads: Record<string, string[]> = {};
    const atoms: AtomProxy[] = [];
    let atomname: string;
    let index: number;
    for (const bead of collection.beads) {
        output += bead.name + " ";
        for (const atom of bead.atoms) {
            atomname = collection.atomName(atom);
            if (atomToBeads[atomname] === undefined) {
                atomToBeads[atomname] = [];
                atoms.push(atom);
            }
            atomToBeads[atomname].push(bead.name ?? "");
        }
    }
    output += "\n\n";

    output += "[ atoms ]\n";
    index = 0;
    atoms.sort(function(a, b) {return a.index - b.index});
    for (const atom of atoms) {
        index += 1;
        atomname = collection.atomName(atom);
        output += index + "\t" + atomname;
        for (const bead of atomToBeads[atomname]) {
            output += "\t" + bead;
        }
        output += "\n";
    }
    return output;
}

/**
 * Shaker (https://github.com/Lp0lp/shaker) Python assignments dict literal:
 * `mapping = { resname: { beadName: {type, charge, atoms: [...]}, ... } }`.
 * The resname used is read from the first bead that actually has atoms
 * assigned (falls back to "UNK" if none do). Atom lists are expanded per
 * bead.atomWeights, so an atom weighted ×2 (pulled toward by repeated
 * clicks) appears twice. This is the exact format parseShakerMapping reads
 * back in for the "Load mapping" feature, so the two must stay in sync.
 * @param collection - BeadCollection
 * @returns empty string if the collection has no beads
 */
export function generatePythonAssignments(collection: BeadCollection): string {
    const beads = collection.beads || [];
    if (beads.length === 0) return "";

    let resname = "UNK";
    for (const bead of beads) {
        if (bead.atoms && bead.atoms.length > 0) { resname = bead.resname; break; }
    }

    const lines: string[] = [];
    lines.push("mapping = {");
    lines.push("    ## resname");
    lines.push(`    "${resname}": {`);
    lines.push("    ## bead name; type(opt.);   charge(opt.);         Mapping.");

    for (const bead of beads) {
        const beadName   = bead.name;
        const beadType   = bead.type || "type";
        const beadCharge = bead.charge ?? 0;
        const atomNames  = collection.expandedAtomNames(bead)
            .map((name) => `'${name}'`).join(",");
        lines.push(
            `        "${beadName}": {"type": "${beadType}", "charge": ${beadCharge}, "atoms": [${atomNames}]},`
        );
    }

    lines.push("    },");
    lines.push("}");
    return lines.join("\n") + "\n";
}

/**
 * PyCGTOOL mapping text: a `[ resname ]` section header followed by one line
 * per bead: `<name> <type> <charge> <atom1> <atom2> ...`. Charge is written as
 * an integer (PyCGTOOL convention). Atoms weighted ×N appear N times —
 * PyCGTOOL's geometric weighting sums duplicate entries, giving them
 * proportionally more influence on the bead centre.
 * @param collection - BeadCollection
 */
export function generatePyCGTOOL(collection: BeadCollection): string {
    const beads = collection.beads;
    if (beads.length === 0) return '';

    let resname = 'UNK';
    for (const bead of beads) {
        if (bead.atoms && bead.atoms.length > 0) { resname = bead.resname; break; }
    }

    const nameW  = Math.max(...beads.map((b) => b.name.length));
    const typeW  = Math.max(...beads.map((b) => (b.type || 'C1').length));
    const chargeW = Math.max(...beads.map((b) => String(Math.round(b.charge ?? 0)).length));

    const lines = [`; CGBuilder export`, `[ ${resname} ]`];
    for (const bead of beads) {
        const atoms  = collection.expandedAtomNames(bead).join(' ');
        const charge = Math.round(bead.charge ?? 0);
        lines.push(
            `${bead.name.padEnd(nameW)}  ${(bead.type || 'C1').padEnd(typeW)}  ${String(charge).padStart(chargeW)}  ${atoms}`
        );
    }
    return lines.join('\n') + '\n';
}

/**
 * Bartender mapping text: one line per bead, `<beadNumber> <idx1>,<idx2>,...`.
 * Bead numbers are 1-based. Atom indices are 1-based (GROMACS convention).
 * Atoms weighted ×N appear N times, matching the Shaker convention.
 * @param collection - BeadCollection
 */
export function generateBartender(collection: BeadCollection): string {
    const lines = ['BEADS'];
    collection.beads.forEach((bead, i) => {
        const indices = bead.expandedAtoms().map((a) => a.index + 1);
        lines.push(`${i + 1} ${indices.join(',')}`);
    });
    return lines.join('\n') + '\n';
}

/**
 * GROMACS .gro text: fixed-column format (resid/resname/atomname/atomid
 * each 5 characters, coordinates 8-character fields), one ATOM line per
 * bead at its geometric centre, coordinates converted from the app's
 * internal Angstrom to .gro's nm (divide by 10).
 *
 * The box vector line is always written as a fixed placeholder ("10 10
 * 10"), not the real extent of the structure — this is a real
 * simplification, not a computed value, since the mapping tool has no
 * concept of a simulation box. Callers piping this into GROMACS should
 * replace that line with an appropriate box for their system.
 * @param collection - BeadCollection
 */
export function generateGRO(collection: BeadCollection): string {
    let output = "Generated with cgbuilder\n" + collection.beads.length + "\n";
    let counter = 0;
    for (const bead of collection.beads) {
        counter += 1;
        const resid    = String(bead.resid).padStart(5);
        const atomid   = String(counter).padStart(5);
        const resname  = bead.resname.padEnd(5).substring(0, 5);
        const atomname = (bead.name ?? "").padStart(5).substring(0, 5);
        const center   = bead.center;
        const x = (center.x / 10).toFixed(3).padStart(8);
        const y = (center.y / 10).toFixed(3).padStart(8);
        const z = (center.z / 10).toFixed(3).padStart(8);
        output += resid + resname + atomname + atomid + x + y + z + '\n';
    }
    output += "10 10 10\n";
    return output;
}

/**
 * Trigger a browser file download of plain text, with no server involved —
 * a throwaway `<a download>` link clicked programmatically. Adapted from
 * <https://ourcodeworld.com/articles/read/189/how-to-create-a-file-and-generate-a-download-with-javascript-in-the-browser-without-a-server>.
 * @param filename
 * @param text
 */
export function download(filename: string, text: string): void {
    const element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

/**
 * Copy text to the clipboard, using the async Clipboard API when available
 * (requires a secure context — HTTPS or localhost) and falling back to the
 * legacy `document.execCommand("copy")` trick (a hidden, focused, selected
 * textarea) otherwise, e.g. when the app is served over plain HTTP.
 * @param text
 */
export function copyTextToClipboard(text: string): Promise<void> {
    if (navigator.clipboard && window.isSecureContext) {
        return navigator.clipboard.writeText(text);
    }
    const textarea = document.createElement("textarea");
    textarea.value = text;
    textarea.style.position = "fixed";
    textarea.style.left = "-9999px";
    textarea.style.top = "-9999px";
    document.body.appendChild(textarea);
    textarea.focus();
    textarea.select();
    try {
        document.execCommand("copy");
    } finally {
        document.body.removeChild(textarea);
    }
    return Promise.resolve();
}

/**
 * Atom names from a PDB/PQR-style file, read directly from columns 13-16
 * of every ATOM/HETATM record, in file order. Bypasses NGL's own structure
 * parsing entirely — this exists only to capture the EXACT original name
 * strings (NGL's parsed `atom.atomname` isn't guaranteed to match the
 * source file byte-for-byte), for BeadCollection.atomName's index-keyed
 * lookup (see setOriginalAtomNames in bead.ts) to fall back to.
 * @param content - raw file text
 * @returns one name per ATOM/HETATM line, in file order
 */
function parsePDBAtomNames(content: string): string[] {
    const names: string[] = [];
    const lines = content.split(/\r?\n/);
    for (const line of lines) {
        if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
            names.push(line.substring(12, 16).trim());
        }
    }
    return names;
}

/**
 * Atom names from a .gro file, read directly from columns 11-15 of each
 * atom line (same rationale as parsePDBAtomNames — capture the original
 * strings independent of NGL's own parsing). The atom count on line 2
 * bounds how many of the following lines are read as atom records.
 * @param content - raw file text
 * @returns one name per atom line, in file order
 */
function parseGROAtomNames(content: string): string[] {
    const names: string[] = [];
    const lines = content.split(/\r?\n/);
    if (lines.length < 3) return names;
    const count = parseInt(lines[1].trim(), 10);
    if (!Number.isFinite(count) || count <= 0) return names;
    const atomLines = lines.slice(2, 2 + count);
    for (const line of atomLines) {
        if (line.length >= 15) names.push(line.substring(10, 15).trim());
    }
    return names;
}

/**
 * Dispatch to the right original-atom-name parser by file extension.
 * @param content - raw file text
 * @param filename - used only for its extension
 * @returns empty array for unrecognised extensions
 */
function parseOriginalAtomNames(content: string, filename: string): string[] {
    const lower = (filename || "").toLowerCase();
    if (lower.endsWith(".pdb") || lower.endsWith(".ent")) return parsePDBAtomNames(content);
    if (lower.endsWith(".gro")) return parseGROAtomNames(content);
    return [];
}

/**
 * Read a just-uploaded File's original atom names (see
 * parseOriginalAtomNames), for BeadCollection.setOriginalAtomNames.
 * Gzipped files are skipped (returns no names, not an error) since this
 * reads plain text directly rather than decompressing. Never rejects —
 * any read/parse failure resolves to an empty array, since this is a
 * best-effort enhancement (falling back to NGL's own atom names) rather
 * than something the rest of loading should fail on.
 * @param file
 */
export function readOriginalAtomNames(file: File): Promise<string[]> {
    const lower = (file && file.name ? file.name : "").toLowerCase();
    if (lower.endsWith(".gz")) return Promise.resolve([]);
    return file.text()
        .then((content) => parseOriginalAtomNames(content, file.name))
        .catch(() => []);
}

/**
 * Parse a pasted Shaker assignments dict literal (the format
 * generatePythonAssignments produces, see its comment) back into bead
 * definitions, for the "Load mapping" feature. Regex-based rather than a
 * real Python/JSON parser — deliberately tolerant of the exact dict-literal
 * syntax Shaker uses (single-quoted atom names, a bare resname wrapper)
 * rather than requiring strict JSON. Bead order in the result follows
 * occurrence order in the text, not any field in the dict itself.
 *
 * A bead's own `{...}` body is matched first (deliberately excluding `{`/`}`
 * from its contents — the thing that distinguishes ONE bead's flat
 * type/charge/atoms object from the surrounding resname wrapper, which
 * nests bead entries inside ITS OWN braces and so can never match this
 * pattern itself), then type/charge/atoms are each pulled out of that body
 * independently. A Python dict literal's key order has no meaning to
 * Shaker, so this doesn't assume "type" precedes "charge" precedes "atoms"
 * the way a single combined pattern would.
 * @param text - pasted Shaker `mapping = {...}` text
 */
export function parseShakerMapping(text: string): BeadDef[] {
    const beads: BeadDef[] = [];
    const beadRe = /"([^"]+)":\s*\{([^{}]*)\}/g;
    let m: RegExpExecArray | null;
    while ((m = beadRe.exec(text)) !== null) {
        const body = m[2];
        const typeMatch = body.match(/"type":\s*"([^"]*)"/);
        const chargeMatch = body.match(/"charge":\s*([^,\s}]+)/);
        const atomsMatch = body.match(/"atoms":\s*\[([^\]]*)\]/);
        if (!typeMatch || !chargeMatch || !atomsMatch) continue;
        const atoms = [...atomsMatch[1].matchAll(/['"]([^'"]+)['"]/g)].map(x => x[1]);
        beads.push({ name: m[1], type: typeMatch[1], charge: parseFloat(chargeMatch[1]) || 0, atoms });
    }
    return beads;
}

/**
 * NGL "ball+stick" representation params with multiple-bond rendering
 * turned on, so a resolved double/triple bond (see chemistry.ts's
 * perceiveChemistry, written back onto the structure's BondProxy.bondOrder)
 * actually displays as parallel sticks instead of NGL's default single
 * line regardless of order. `overrides` is merged in on top, last, so
 * callers can adjust individual params (e.g. color, opacity) without
 * losing the multiple-bond defaults.
 * @param overrides
 * @returns NGL representation parameters
 */
export function bondAwareRepresentationParams(overrides: RepresentationParams = {}): RepresentationParams {
    return Object.assign({ multipleBond: true, bondSpacing: 1, bondScale: 0.4 }, overrides);
}
