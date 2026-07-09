/* ===========================================================================
   RDKit.js loader
   ===========================================================================
   Loads RDKit_minimal lazily on first use. Bead-type prediction is the only
   feature that needs it, so most page loads never pay for it at all. */

let _rdkitPromise = null;

/**
 * Load RDKit_minimal from the unpkg CDN and initialize its WASM module,
 * caching the result so subsequent calls return the same resolved instance
 * instantly rather than re-injecting the script tag. Also reuses
 * `window.RDKit` if some other code already loaded it first. A failed
 * attempt (network error, bad WASM init) clears the cache rather than
 * caching the rejection, so the next call retries from scratch instead of
 * re-rejecting the same dead promise for the rest of the page session.
 * @returns {Promise<object>} resolves to the initialized RDKit module
 */
export function loadRDKit() {
    if (_rdkitPromise) return _rdkitPromise;
    _rdkitPromise = new Promise((resolve, reject) => {
        if (window.RDKit) { resolve(window.RDKit); return; }
        const script = document.createElement('script');
        script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
        script.onload = () => {
            window.initRDKitModule()
                .then(rdkit => { window.RDKit = rdkit; resolve(rdkit); })
                .catch(reject);
        };
        script.onerror = () => reject(new Error('Failed to load RDKit.js'));
        document.head.appendChild(script);
    }).catch((err) => {
        _rdkitPromise = null;
        throw err;
    });
    return _rdkitPromise;
}
