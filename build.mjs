/* ===========================================================================
   esbuild build driver
   ===========================================================================
   Bundles src/main.ts (and its imports — including the NGL npm package) into
   a single ESM file at dist/main.js, then copies the static public/ assets
   over it. Type-checking is a separate step (`tsc --noEmit`, run by the
   `build` npm script before this); esbuild only strips types and bundles.
   Note: RDKit and marked are still loaded as browser globals (RDKit via a
   runtime <script> injection in rdkit.ts, marked via a CDN <script> tag in
   index.html), so they are intentionally NOT bundled here.

   With --watch, also starts an HTTP dev server on dist/ with live reload:
   a rebuild re-copies public/ and pushes a browser reload via SSE. */
import * as esbuild from "esbuild";
import { cpSync, rmSync } from "node:fs";

const watch = process.argv.includes("--watch");

/** Copy public's contents into dist, alongside the freshly-built bundle.
    Mirrors the old `cp -r public/. dist/`. Runs after every (re)build so
    static assets stay in sync in watch mode too. */
const copyPublicPlugin = {
  name: "copy-public",
  setup(build) {
    build.onEnd(() => {
      cpSync("public", "dist", { recursive: true });
    });
  },
};

/** @type {import('esbuild').BuildOptions} */
const options = {
  entryPoints: ["src/main.ts"],
  bundle: true,
  format: "esm",
  target: "es2022",
  platform: "browser",
  sourcemap: true,
  outdir: "dist",
  // Code-splitting: the dynamic import('@rdkit/rdkit') in rdkit.ts is emitted
  // as its own chunk, loaded only on first bead-type prediction. Requires
  // outdir (not outfile). The entry keeps its name -> dist/main.js.
  splitting: true,
  // RDKit's WASM binary is copied out as a build asset; the import yields its
  // (page-relative) URL, handed to RDKit via locateFile (see rdkit.ts).
  loader: { ".wasm": "file" },
  // RDKit's Emscripten glue references Node built-ins (fs/crypto/...) inside
  // dead `if (ENVIRONMENT_IS_NODE)` branches that never run in the browser.
  // Mark them external so esbuild leaves those unreachable require() calls
  // alone instead of trying (and failing) to resolve them.
  external: ["fs", "path", "crypto", "module", "worker_threads"],
  logLevel: "info",
  plugins: [copyPublicPlugin],
  // Only in dev: subscribe to esbuild's SSE endpoint and reload on rebuild.
  ...(watch && {
    banner: {
      js: `new EventSource('/esbuild').addEventListener('change', () => location.reload());`,
    },
  }),
};

rmSync("dist", { recursive: true, force: true });

if (watch) {
  const ctx = await esbuild.context(options);
  await ctx.watch(); // triggers a rebuild (and copyPublic via onEnd) on changes
  const { host, port } = await ctx.serve({ servedir: "dist" });
  console.log(
    `esbuild: serving http://localhost:${port} (watching for changes…)`,
  );
} else {
  await esbuild.build(options); // onEnd copies public
}
