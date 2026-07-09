/* ===========================================================================
   Ambient browser globals
   ===========================================================================
   NGL, RDKit and marked are all loaded via plain <script> tags in
   index.html rather than ES imports, so they exist as browser globals with
   no module of their own to import types from. Declare just the surface this
   app touches. */
import type { NGLStatic, RDKitModule } from './types.js';

declare global {
    /** NGL viewer library (https://nglviewer.org). */
    const NGL: NGLStatic;

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
