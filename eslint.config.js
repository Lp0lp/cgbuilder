const js = require('@eslint/js');

module.exports = [
    {
        ignores: ['eslint.config.js', 'prettier.config.js'],
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
                // External library loaded via CDN
                NGL: 'readonly',
            },
        },
    },
];
