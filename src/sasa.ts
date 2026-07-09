/* ===========================================================================
   Solvent-accessible surface area (SASA)
   ===========================================================================
   Numerical SASA via the Shrake-Rupley algorithm (Shrake & Rupley, 1973,
   https://doi.org/10.1016/0022-2836(73)90011-9 ): for every atom/bead, sample
   points on a sphere expanded by the solvent probe radius, count the
   fraction not buried inside any neighbour's own probe-expanded sphere, and
   scale by that sphere's area. Summed over all particles, this gives a
   numerical estimate of the true analytical SASA without needing exact
   sphere-sphere intersection geometry (see shrakeRupley for the per-atom
   detail).

   This is computed independently for the all-atom structure (aaSASA) and
   for the coarse-grained bead representation (cgSASA), using the SAME probe
   radius and algorithm for both, so the resulting Å² values are directly
   comparable — that comparison (the Mapping panel's SASA Δ) is the actual
   point of this module. That numerical computation is entirely independent
   of NGL's own "surface" representation (the translucent mesh toggled in
   the viewer) — every Å² value comes from the hand-implemented shrakeRupley
   below, not from NGL.

   This file also bridges to that NGL visual surface, though, via
   beadsToPDB: it reuses the same per-bead radius (BEAD_RADII/beadRadius)
   the numerical cgSASA uses, but writes it into a synthetic PDB's
   temperature-factor column so NGL's renderer can read individual bead
   radii back out (see beadsToPDB's own comment for the mechanism). So the
   two halves of this file share their radius data, even though the
   numbers and the visual mesh are computed by entirely different code.

   Probe radius in Angstrom. 0.191 nm = 1.91 Å is the Martini tiny-bead radius,
   used here so the AA and CG surfaces are computed with the same probe and are
   directly comparable. */
import type { Bead, BeadCollection } from './bead.js';
import type { Structure } from './types.js';

/** A particle for shrakeRupley: [x, y, z, vdwRadius] in Angstrom. */
type Particle = [number, number, number, number];

export const PROBE_RADIUS = 1.91;

/* Martini bead vdW radii in Angstrom, keyed by size class. Values come from the
   standard Martini bead radii (nm -> A): regular 0.264, small 0.230, tiny 0.191.
   "U" beads are virtual/ghost beads and contribute no surface (radius 0).
   Edit this table to switch force fields. */
export const BEAD_RADII: Record<string, number> = {
    R: 2.64,
    S: 2.30,
    T: 1.91,
    U: 0.0,
};
const DEFAULT_BEAD_SIZE = "R";

/**
 * The Martini size-class letter (R/S/T/U) implied by a bead type string,
 * read from its first character — e.g. "SP2a" -> "S", "TC5" -> "T", a plain
 * "P2" (no size prefix) -> the "R" (regular) default.
 * @param type - bead type code, e.g. "SP2a"
 * @returns one of "R", "S", "T", "U"
 */
export function beadSizeClass(type: string): string {
    if (!type) return DEFAULT_BEAD_SIZE;
    const first = type.trim().charAt(0).toUpperCase();
    if (first === "S" || first === "T" || first === "U") return first;
    return DEFAULT_BEAD_SIZE;
}

/**
 * The Martini vdW radius (Å) for a bead, from its type's size class.
 * @param bead - has a `type` string (e.g. "SP2a")
 * @returns radius in Angstrom (0 for virtual "U" beads)
 */
export function beadRadius(bead: Bead): number {
    return BEAD_RADII[beadSizeClass(bead.type)];
}

// Standard vdW radii (Angstrom) for elements found in organic molecules.
const VDW_RADII: Record<string, number> = {
    H: 1.09, C: 1.75, N: 1.61, O: 1.56, F: 1.44,
    P: 1.80, S: 1.79, CL: 1.74, BR: 1.85, I: 2.00,
};
const DEFAULT_VDW_RADIUS = 1.75;

/**
 * Fibonacci golden-ratio lattice — n near-uniformly distributed points on
 * the unit sphere, used by shrakeRupley as the fixed set of sample
 * directions tested around every atom. Deterministic (no RNG), so repeated
 * calls with the same n always sample the same directions.
 * @param n - number of points to generate
 * @returns unit-length [x,y,z] points
 */
function fibonacciSpherePoints(n: number): [number, number, number][] {
    const pts: [number, number, number][] = new Array(n);
    const phi = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < n; i++) {
        const y = 1 - (i / (n - 1)) * 2;
        const r = Math.sqrt(1 - y * y);
        const theta = phi * i;
        pts[i] = [r * Math.cos(theta), y, r * Math.sin(theta)];
    }
    return pts;
}

/**
 * Numerical SASA (Å²) via Shrake-Rupley: for each particle i, place nPoints
 * test points on a sphere of radius (vdwRadius_i + probeRadius) centred on
 * it — this is the surface the centre of a solvent probe sphere could touch
 * while still contacting atom i. A test point is "buried" if it falls
 * inside ANY other particle j's own probe-expanded sphere (i.e. the probe
 * centred there would simultaneously overlap j), meaning solvent can't
 * actually reach that point while resting against i. The fraction of test
 * points that are NOT buried, times the full sphere area, is i's
 * contribution; summing over every particle gives the total SASA.
 *
 * The inner neighbour search first prunes to particles whose probe-expanded
 * spheres could possibly overlap i's own (a generous distance cutoff),
 * before testing individual sample points against only that shortlist —
 * an O(n²) neighbour pass plus O(n_neighbours × nPoints) point tests per
 * particle, not a full O(n² × nPoints).
 * @param particles - [x,y,z,vdwRadius] in Angstrom
 * @param probeRadius - solvent probe radius in Angstrom
 * @param nPoints - sample points per particle (more = smoother, slower)
 * @returns total SASA in Å²
 */
export function shrakeRupley(particles: Particle[], probeRadius: number, nPoints = 4800): number {
    const n = particles.length;
    if (n === 0) return 0;
    const unitPts = fibonacciSpherePoints(nPoints);
    let totalSASA = 0;

    for (let i = 0; i < n; i++) {
        const [xi, yi, zi, ri] = particles[i];
        const shellR = ri + probeRadius;

        // Each neighbour's coordinates/cutoff are fixed regardless of which
        // sample point is being tested, so they're resolved once here rather
        // than re-derived from `particles[j]` on every one of the nPoints
        // iterations below.
        const neighbors: Particle[] = [];
        for (let j = 0; j < n; j++) {
            if (j === i) continue;
            const [xj, yj, zj, rj] = particles[j];
            const cutoff = shellR + rj + probeRadius;
            const dx = xi - xj, dy = yi - yj, dz = zi - zj;
            if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) neighbors.push([xj, yj, zj, rj + probeRadius]);
        }

        let exposed = 0;
        for (const [ux, uy, uz] of unitPts) {
            const px = xi + shellR * ux;
            const py = yi + shellR * uy;
            const pz = zi + shellR * uz;
            let buried = false;
            for (const [xj, yj, zj, cutoff] of neighbors) {
                const dx = px - xj, dy = py - yj, dz = pz - zj;
                if (dx*dx + dy*dy + dz*dz < cutoff*cutoff) { buried = true; break; }
            }
            if (!buried) exposed++;
        }
        totalSASA += (exposed / nPoints) * 4 * Math.PI * shellR * shellR;
    }
    return totalSASA;
}

/**
 * Whole-structure SASA (Å²) for the all-atom representation, using
 * vdW radii per element.
 * @param structure - NGL-style structure (eachAtom)
 * @param probeRadius - solvent probe radius in Angstrom
 */
export function aaSASA(structure: Structure, probeRadius: number): number {
    const particles: Particle[] = [];
    structure.eachAtom((atom) => {
        const el = (atom.element || "").toUpperCase();
        const r = VDW_RADII[el] ?? DEFAULT_VDW_RADIUS;
        particles.push([atom.x, atom.y, atom.z, r]);
    });
    return shrakeRupley(particles, probeRadius);
}

/**
 * Whole-mapping SASA (Å²) for the coarse-grained bead representation, using
 * each bead's own Martini size-class radius (see beadRadius) centred at its
 * geometric centre. Empty beads and zero-radius "U" (virtual) beads
 * contribute no surface and are skipped.
 * @param collection - BeadCollection (has a `beads` array)
 * @param probeRadius - solvent probe radius in Angstrom (same value as used
 *   for aaSASA, so the two totals are comparable)
 */
export function cgSASA(collection: BeadCollection, probeRadius: number): number {
    const particles: Particle[] = [];
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) continue;
        const r = beadRadius(bead);
        if (r <= 0) continue;
        const c = bead.center;
        particles.push([c.x, c.y, c.z, r]);
    }
    return shrakeRupley(particles, probeRadius);
}

/**
 * Pad an atom name into PDB columns 13-16.
 * @param name
 * @returns exactly 4 characters
 */
function formatPDBAtomName(name: string): string {
    name = (name || "").substring(0, 4);
    if (name.length >= 4) return name;
    return (" " + name).padEnd(4);
}

/**
 * Serialise the CG beads to a minimal PDB string, one ATOM per bead, with
 * the bead's own Martini vdW radius stored in the temperature-factor column
 * (cols 61-66). Exists so the CG representation can be loaded into NGL as a
 * real structure and given its own "surface" representation (the visual
 * toggle in the viewer) the same way the all-atom structure has — NGL has
 * no native concept of a "bead", so this fakes a PDB it already knows how
 * to render.
 *
 * The temperature-factor placement specifically is the trick that makes
 * per-bead radii (R/S/T are different sizes) work in NGL's surface
 * renderer at all: NGL's "surface" representation normally derives radius
 * from the atom's element, which would force every bead to the same size.
 * Passing `radiusType: "bfactor"` in the representation params (see
 * drawCGSurface in visualization.ts, the only caller of this function)
 * tells NGL to read each atom's per-bead radius straight back out of the
 * column this function just wrote it into, instead.
 *
 * Beads with no atoms or a zero radius (virtual "U" beads) are skipped.
 * @param collection - BeadCollection (has a `beads` array)
 * @returns PDB-format text (empty string if no beads qualify)
 */
export function beadsToPDB(collection: BeadCollection): string {
    const lines: string[] = [];
    let serial = 0;
    for (const bead of collection.beads) {
        if (bead.atoms.length === 0) continue;
        const radius = beadRadius(bead);
        if (radius <= 0) continue;
        serial += 1;
        const center = bead.center;
        const serStr = String(serial % 100000).padStart(5);
        const name = formatPDBAtomName(bead.name ?? "");
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
