/* ===========================================================================
   Ambient browser globals
   ===========================================================================
   RDKit is still loaded via a plain <script> tag injected at runtime by
   rdkit.ts rather than an ES import, so it exists as a browser global with no
   module of its own to import types from. Declare just the surface this app
   touches.

   (NGL and marked used to be here too; both are now npm imports — NGL adapted
   in ngl.ts, marked imported directly in main.ts.) */
import type { RDKitModule } from './types.js';

declare global {
    interface Window {
        /** Cached RDKit_minimal module once loaded (see rdkit.ts). */
        RDKit?: RDKitModule;
        /** RDKit_minimal's WASM initializer, injected by its script tag. */
        initRDKitModule(): Promise<RDKitModule>;
    }
}

export {};
