const js = require('@eslint/js');

module.exports = [
    {
        ignores: ['eslint.config.js', 'prettier.config.js', 'vitest.config.js'],
    },
    js.configs.recommended,
    {
        files: ['scripts/**/*.js'],
        languageOptions: {
            ecmaVersion: 2020,
            sourceType: 'script',
            globals: {
                // Browser globals
                window: 'readonly',
                document: 'readonly',
                navigator: 'readonly',
                console: 'readonly',
                setTimeout: 'readonly',
                Promise: 'readonly',
                module: 'readonly',
                // External library loaded via CDN
                NGL: 'readonly',
            },
        },
    },
    {
        files: ['tests/**/*.js'],
        languageOptions: {
            ecmaVersion: 2020,
            sourceType: 'script',
            globals: {
                describe: 'readonly',
                it: 'readonly',
                expect: 'readonly',
                beforeEach: 'readonly',
                afterEach: 'readonly',
                vi: 'readonly',
                require: 'readonly',
                __dirname: 'readonly',
            },
        },
    },
];
