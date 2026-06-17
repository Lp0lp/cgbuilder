let _rdkitPromise = null;

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
    });
    return _rdkitPromise;
}
