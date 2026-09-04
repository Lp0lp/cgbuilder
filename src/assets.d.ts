/* ===========================================================================
   Asset module declarations
   ===========================================================================
   Lets TypeScript understand esbuild's `file` loader imports, where importing
   a binary asset yields the URL string esbuild emits for it. Used by rdkit.ts
   to locate the RDKit WASM binary at runtime (see its `locateFile`). */
declare module '*.wasm' {
  const url: string;
  export default url;
}
