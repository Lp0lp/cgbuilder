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
            if (!atomToBeads[atomname]) {
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

/* Taken from <https://ourcodeworld.com/articles/read/189/how-to-create-a-file-and-generate-a-download-with-javascript-in-the-browser-without-a-server> */
export function download(filename, text) {
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

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

function parseOriginalAtomNames(content, filename) {
    const lower = (filename || "").toLowerCase();
    if (lower.endsWith(".pdb") || lower.endsWith(".ent")) return parsePDBAtomNames(content);
    if (lower.endsWith(".gro")) return parseGROAtomNames(content);
    return [];
}

export function readOriginalAtomNames(file) {
    const lower = (file && file.name ? file.name : "").toLowerCase();
    if (lower.endsWith(".gz")) return Promise.resolve([]);
    return file.text()
        .then((content) => parseOriginalAtomNames(content, file.name))
        .catch(() => []);
}

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

export function bondAwareRepresentationParams(overrides = {}) {
    return Object.assign({ multipleBond: true, bondSpacing: 1, bondScale: 0.4 }, overrides);
}
