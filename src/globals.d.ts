/* ===========================================================================
   Ambient browser globals
   ===========================================================================
   RDKit and marked are still loaded via plain <script> tags (RDKit injected
   at runtime by rdkit.ts, marked via a CDN tag in index.html) rather than ES
   imports, so they exist as browser globals with no module of their own to
   import types from. Declare just the surface this app touches.

   (NGL is no longer here — it's now the `ngl` npm package, adapted in
   ngl.ts.) */
import type { RDKitModule } from './types.js';

declare global {
    /** marked markdown renderer, used for the doc pages (see main.ts). */
    const marked: { parse(markdown: string): string };

    interface Window {
        /** Cached RDKit_minimal module once loaded (see rdkit.ts). */
        RDKit?: RDKitModule;
        /** RDKit_minimal's WASM initializer, injected by its script tag. */
        initRDKitModule(): Promise<RDKitModule>;
    }
}

export {};
