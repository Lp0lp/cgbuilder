/* ===========================================================================
   esbuild build driver
   ===========================================================================
   Bundles src/main.ts (and its imports — including the NGL npm package) into
   a single ESM file at dist/main.js, then copies the static public/ assets
   over it. Type-checking is a separate step (`tsc --noEmit`, run by the
   `build` npm script before this); esbuild only strips types and bundles.

   Note: RDKit and marked are still loaded as browser globals (RDKit via a
   runtime <script> injection in rdkit.ts, marked via a CDN <script> tag in
   index.html), so they are intentionally NOT bundled here. */
import * as esbuild from 'esbuild';
import { cpSync, rmSync } from 'node:fs';

const watch = process.argv.includes('--watch');

/** @type {import('esbuild').BuildOptions} */
const options = {
    entryPoints: ['src/main.ts'],
    bundle: true,
    format: 'esm',
    target: 'es2022',
    platform: 'browser',
    sourcemap: true,
    outfile: 'dist/main.js',
    logLevel: 'info',
};

function copyPublic() {
    // Mirror the old `cp -r public/. dist/`: copy public's contents into dist,
    // alongside the freshly-built bundle.
    cpSync('public', 'dist', { recursive: true });
}

rmSync('dist', { recursive: true, force: true });

if (watch) {
    const ctx = await esbuild.context(options);
    await ctx.rebuild();
    copyPublic();
    await ctx.watch();
    console.log('esbuild: watching for changes…');
} else {
    await esbuild.build(options);
    copyPublic();
}
