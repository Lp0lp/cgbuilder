import { defineConfig } from 'vitest/config';

// Scope the test run to the current TypeScript suite under src/tests only.
// The archived plain-JS snapshot in previous/ ships its own duplicate copy of
// these tests (previous/tests/*.test.js); without this, `vitest` would pick
// those up too and run the same suite twice.
export default defineConfig({
    test: {
        include: ['src/tests/**/*.test.ts'],
    },
});
