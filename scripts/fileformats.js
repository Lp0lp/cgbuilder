/* ===========================================================================
   File I/O — export generators, import parsers, small browser helpers
   ===========================================================================
   Three unrelated concerns share this file because each is small enough not
   to need its own module:
     1. Output generators (generateNDX, generateMap, generatePythonAssignments,
        generateGRO) — serialise the current bead mapping into the formats
        downstream tools expect (GROMACS .ndx/.gro, a Martini .map, a Shaker
        Python assignments dict).
     2. Import parsers (parsePDBAtomNames/parseGROAtomNames/
        parseOriginalAtomNames/readOriginalAtomNames, parseShakerMapping) —
        the reverse direction: recovering the original AA atom names from a
        just-loaded file (NGL's own parsed atom names aren't guaranteed to
        match the source file byte-for-byte), and re-importing a previously
        exported Shaker mapping.
     3. Small browser utility helpers (download, copyTextToClipboard) and one
        NGL representation-params helper (bondAwareRepresentationParams) that
        don't have a more specific home elsewhere. */

/**
 * GROMACS .ndx text: one `[ bead name ]` group per bead, listing its atoms'
 * 1-based (GROMACS convention) indices.
 * @param {object} collection - BeadCollection
 * @returns {string}
 */
export function generateNDX(collection) {
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
 * @param {object} collection - BeadCollection
 * @returns {string}
 */
export function generateMap(collection) {
    let output = "[ to ]\nmartini\n\n[ martini ]\n";
    let atomToBeads = {};
    let atoms = [];
    let atomname;
    let index;
    for (const bead of collection.beads) {
        output += bead.name + " ";
        for (const atom of bead.atoms) {
            atomname = collection.atomName(atom);
            if (atomToBeads[atomname] === undefined) {
                atomToBeads[atomname] = [];
                atoms.push(atom);
            }
            atomToBeads[atomname].push(bead.name);
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
 * @param {object} collection - BeadCollection
 * @returns {string} empty string if the collection has no beads
 */
export function generatePythonAssignments(collection) {
    const beads = collection.beads || [];
    if (beads.length === 0) return "";

    let resname = "UNK";
    for (const bead of beads) {
        if (bead.atoms && bead.atoms.length > 0) { resname = bead.resname; break; }
    }

    let lines = [];
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
 * @param {object} collection - BeadCollection
 * @returns {string}
 */
export function generateGRO(collection) {
    let resid = "    0";
    let resname = "";
    let atomname = "    0";
    let atomid = 0;
    let x, y, z, center;
    let output = "Generated with cgbuilder\n" + collection.beads.length + "\n";
    let counter = 0;
    for (const bead of collection.beads) {
        counter += 1;
        resid    = String(bead.resid).padStart(5);
        atomid   = String(counter).padStart(5);
        resname  = bead.resname.padEnd(5).substring(0, 5);
        atomname = bead.name.padStart(5).substring(0, 5);
        center   = bead.center;
        x = (center.x / 10).toFixed(3).padStart(8);
        y = (center.y / 10).toFixed(3).padStart(8);
        z = (center.z / 10).toFixed(3).padStart(8);
        output += resid + resname + atomname + atomid + x + y + z + '\n';
    }
    output += "10 10 10\n";
    return output;
}

/**
 * Trigger a browser file download of plain text, with no server involved —
 * a throwaway `<a download>` link clicked programmatically. Adapted from
 * <https://ourcodeworld.com/articles/read/189/how-to-create-a-file-and-generate-a-download-with-javascript-in-the-browser-without-a-server>.
 * @param {string} filename
 * @param {string} text
 */
export function download(filename, text) {
    let element = document.createElement('a');
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
 * @param {string} text
 * @returns {Promise<void>}
 */
export function copyTextToClipboard(text) {
    if (navigator.clipboard && window.isSecureContext) {
        return navigator.clipboard.writeText(text);
    }
    let textarea = document.createElement("textarea");
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
 * lookup (see setOriginalAtomNames in bead.js) to fall back to.
 * @param {string} content - raw file text
 * @returns {Array<string>} one name per ATOM/HETATM line, in file order
 */
function parsePDBAtomNames(content) {
    const names = [];
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
 * @param {string} content - raw file text
 * @returns {Array<string>} one name per atom line, in file order
 */
function parseGROAtomNames(content) {
    const names = [];
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
 * @param {string} content - raw file text
 * @param {string} filename - used only for its extension
 * @returns {Array<string>} empty array for unrecognised extensions
 */
function parseOriginalAtomNames(content, filename) {
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
 * @param {File} file
 * @returns {Promise<Array<string>>}
 */
export function readOriginalAtomNames(file) {
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
 * @param {string} text - pasted Shaker `mapping = {...}` text
 * @returns {Array<{name: string, type: string, charge: number, atoms: Array<string>}>}
 */
export function parseShakerMapping(text) {
    const beads = [];
    // Matches each bead line: "NAME": {"type": "T", "charge": 0, "atoms": ['A1', 'A2']}
    const re = /"([^"]+)":\s*\{\s*"type":\s*"([^"]*)",\s*"charge":\s*([^,\s\n]+),\s*"atoms":\s*\[([^\]]*)\]/g;
    let m;
    while ((m = re.exec(text)) !== null) {
        const atoms = [...m[4].matchAll(/['"]([^'"]+)['"]/g)].map(x => x[1]);
        beads.push({ name: m[1], type: m[2], charge: parseFloat(m[3]) || 0, atoms });
    }
    return beads;
}

/**
 * NGL "ball+stick" representation params with multiple-bond rendering
 * turned on, so a resolved double/triple bond (see chemistry.js's
 * perceiveChemistry, written back onto the structure's BondProxy.bondOrder)
 * actually displays as parallel sticks instead of NGL's default single
 * line regardless of order. `overrides` is merged in on top, last, so
 * callers can adjust individual params (e.g. color, opacity) without
 * losing the multiple-bond defaults.
 * @param {object} [overrides]
 * @returns {object} NGL representation parameters
 */
export function bondAwareRepresentationParams(overrides = {}) {
    return Object.assign({ multipleBond: true, bondSpacing: 1, bondScale: 0.4 }, overrides);
}
