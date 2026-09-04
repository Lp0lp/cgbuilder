/* ===========================================================================
   RDKit.js loader
   ===========================================================================
   Loads RDKit_minimal lazily on first use via a dynamic import(), so esbuild
   code-splits the (multi-MB) RDKit module into its own chunk that is only
   fetched when bead-type prediction is actually used — most page loads never
   pay for it at all. The WASM binary is emitted as a build asset (esbuild's
   `file` loader) and its URL is handed to RDKit through `locateFile`. */
import wasmUrl from '@rdkit/rdkit/dist/RDKit_minimal.wasm';
import type { RDKitLoader, RDKitModule } from '@rdkit/rdkit';

let _rdkitPromise: Promise<RDKitModule> | null = null;

/**
 * Dynamically import RDKit_minimal and initialize its WASM module, caching
 * the result so subsequent calls return the same resolved instance instantly
 * rather than re-importing/re-initializing. A failed attempt (chunk load
 * error, bad WASM init) clears the cache rather than caching the rejection,
 * so the next call retries from scratch instead of re-rejecting the same
 * dead promise for the rest of the page session.
 * @returns resolves to the initialized RDKit module
 */
export function loadRDKit(): Promise<RDKitModule> {
  if (_rdkitPromise) return _rdkitPromise;
  _rdkitPromise = import('@rdkit/rdkit')
    // The package is CommonJS (`module.exports = initRDKitModule`), so the
    // loader is the default export at runtime — the shipped .d.ts only
    // declares the named types, hence the cast to reach `.default`.
    .then((mod) => {
      const initRDKitModule = (mod as unknown as { default: RDKitLoader })
        .default;
      return initRDKitModule({ locateFile: () => wasmUrl });
    })
    .catch((err) => {
      _rdkitPromise = null;
      throw err;
    });
  return _rdkitPromise;
}
