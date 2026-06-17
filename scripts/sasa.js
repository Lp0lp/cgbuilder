/* ===========================================================================
   Solvent-accessible surface area (SASA)
   ===========================================================================
   Probe radius in Angstrom. 0.191 nm = 1.91 Å is the Martini tiny-bead radius,
   used here so the AA and CG surfaces are computed with the same probe and are
   directly comparable. */
export const PROBE_RADIUS = 1.91;

/* Martini bead vdW radii in Angstrom, keyed by size class. Values come from the
   standard Martini bead radii (nm -> A): regular 0.264, small 0.230, tiny 0.191.
   "U" beads are virtual/ghost beads and contribute no surface (radius 0).
   Edit this table to switch force fields. */
export const BEAD_RADII = {
    R: 2.64,
    S: 2.30,
    T: 1.91,
    U: 0.0,
};
const DEFAULT_BEAD_SIZE = "R";

export function beadSizeClass(type) {
    if (!type) return DEFAULT_BEAD_SIZE;
    const first = type.trim().charAt(0).toUpperCase();
    if (first === "S" || first === "T" || first === "U") return first;
    return DEFAULT_BEAD_SIZE;
}

export function beadRadius(bead) {
    return BEAD_RADII[beadSizeClass(bead.type)];
}

// Standard Bondi vdW radii (Angstrom) for elements found in organic molecules.
const VDW_RADII = {
    H: 1.09, C: 1.75, N: 1.61, O: 1.56, F: 1.44,
    P: 1.80, S: 1.79, CL: 1.74, BR: 1.85, I: 2.00,
};
const DEFAULT_VDW_RADIUS = 1.75;

// Fibonacci golden-ratio lattice — generates n near-uniformly distributed
// points on the unit sphere.
function fibonacciSpherePoints(n) {
    const pts = new Array(n);
    const phi = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < n; i++) {
        const y = 1 - (i / (n - 1)) * 2;
        const r = Math.sqrt(1 - y * y);
        const theta = phi * i;
        pts[i] = [r * Math.cos(theta), y, r * Math.sin(theta)];
    }
    return pts;
}

// Shrake-Rupley numerical SASA in Å².
// particles: array of [x, y, z, vdwRadius] in Angstrom.
export function shrakeRupley(particles, probeRadius, nPoints = 4800) {
    const n = particles.length;
    if (n === 0) return 0;
    const unitPts = fibonacciSpherePoints(nPoints);
    let totalSASA = 0;

    for (let i = 0; i < n; i++) {
        const [xi, yi, zi, ri] = particles[i];
        const shellR = ri + probeRadius;

        const neighbors = [];
        for (let j = 0; j < n; j++) {
            if (j === i) continue;
            const [xj, yj, zj, rj] = particles[j];
            const cutoff = shellR + rj + probeRadius;
            const dx = xi - xj, dy = yi - yj, dz = zi - zj;
            if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) neighbors.push(j);
        }

        let exposed = 0;
        for (const [ux, uy, uz] of unitPts) {
            const px = xi + shellR * ux;
            const py = yi + shellR * uy;
            const pz = zi + shellR * uz;
            let buried = false;
            for (const j of neighbors) {
                const [xj, yj, zj, rj] = particles[j];
                const cutoff = rj + probeRadius;
                const dx = px - xj, dy = py - yj, dz = pz - zj;
                if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) { buried = true; break; }
            }
            if (!buried) exposed++;
        }
        totalSASA += (exposed / nPoints) * 4 * Math.PI * shellR * shellR;
    }
    return totalSASA;
}

export function aaSASA(structure, probeRadius) {
    const particles = [];
    structure.eachAtom((atom) => {
        const el = (atom.element || "").toUpperCase();
        const r = VDW_RADII[el] ?? DEFAULT_VDW_RADIUS;
        particles.push([atom.x, atom.y, atom.z, r]);
    });
    return shrakeRupley(particles, probeRadius);
}

export function cgSASA(collection, probeRadius) {
    const particles = [];
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) continue;
        const r = beadRadius(bead);
        if (r <= 0) continue;
        const c = bead.center;
        particles.push([c.x, c.y, c.z, r]);
    }
    return shrakeRupley(particles, probeRadius);
}

// Pad an atom name into PDB columns 13-16.
function formatPDBAtomName(name) {
    name = (name || "").substring(0, 4);
    if (name.length >= 4) return name;
    return (" " + name).padEnd(4);
}

// Serialise the CG beads to a minimal PDB string, one ATOM per bead, with the
// bead vdW radius stored in the temperature-factor column (cols 61-66). Beads
// with no atoms or a zero radius (virtual "U" beads) are skipped.
export function beadsToPDB(collection) {
    let lines = [];
    let serial = 0;
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) continue;
        const radius = beadRadius(bead);
        if (radius <= 0) continue;
        serial += 1;
        const center = bead.center;
        const serStr = String(serial % 100000).padStart(5);
        const name = formatPDBAtomName(bead.name);
        const resname = (bead.resname || "BEA").substring(0, 3).padEnd(3);
        const resid = String(((bead.resid % 10000) + 10000) % 10000).padStart(4);
        const x = center.x.toFixed(3).padStart(8);
        const y = center.y.toFixed(3).padStart(8);
        const z = center.z.toFixed(3).padStart(8);
        const occ = "  1.00";
        const bfac = radius.toFixed(2).padStart(6);
        lines.push(
            "ATOM  " + serStr + " " + name + " " + resname + " " + "A" + resid +
            " " + "   " + x + y + z + occ + bfac + "          " + " C"
        );
    }
    if (lines.length === 0) return "";
    lines.push("END");
    return lines.join("\n") + "\n";
}
