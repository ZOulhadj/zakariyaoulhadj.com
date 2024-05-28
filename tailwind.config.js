/** @type {import('tailwindcss').Config} */
export default {
    content: [
        "./src/**/*.{astro,html,js,jsx,md,mdx,svelte,ts,tsx,vue}",
    ],
    theme: {
        extend: {

        },
    },
    // colors: {
    // },
    // typography: (theme) => ({
    //     DEFAULT: {
    //         css: {
    //             color: theme('colors.gray.800'),

    //             // ...
    //         },
    //     },
    // }),
    darkMode: 'class',
    corePlugins: {
        aspectRatio: false,
    },
    plugins: [
        require('@tailwindcss/typography'),
    ],
}
