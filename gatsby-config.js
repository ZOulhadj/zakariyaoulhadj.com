module.exports = {
  siteMetadata: {
    title: `website`,
    description: `Zakariya Oulhadj's personal website`,
    siteUrl: `https://zakariyaoulhadj.com`,
    siteDomain: `zakariyaoulhadj.com`
  },
  plugins: [{
    resolve: "gatsby-plugin-sass",
    options: {
      sassOptions: {
        includePaths: ["node_modules/bootstrap/scss"]
      }
    }
  }, {
    resolve: 'gatsby-plugin-google-gtag',
    options: {
      trackingIds: ["G-8B9M5SYE68"]
    }
  }, "gatsby-plugin-image", "gatsby-plugin-mdx", "gatsby-plugin-sharp", "gatsby-transformer-sharp", {
    resolve: 'gatsby-source-filesystem',
    options: {
      "name": "images",
      "path": "./src/images/"
    },
    __key: "images"
  }, {
    resolve: 'gatsby-source-filesystem',
    options: {
      "name": "pages",
      "path": "./src/pages/"
    },
    __key: "pages"
  }, {
    resolve: 'gatsby-source-filesystem',
    options: {
      "name": "blog",
      "path": "./content/blog/"
    },
    __key: "blog"
  }, {
    resolve: 'gatsby-transformer-remark',
    options: {
      plugins: ['gatsby-remark-autolink-headers']
    }
  }, {
      resolve: 'gatsby-omni-font-loader',
      options: {
        enableListener: true,
        preconnect: ['https://fonts.google.com', 'https://fonts.gstatic.com'],
        web: [
          {
            name: 'FigTree',
            file: 'https://fonts.googleapis.com/css2?family=Figtree:wght@400;700;900&display=swap'
          },
          {
            name: 'Playfair Display',
            file: 'https://fonts.googleapis.com/css2?family=Playfair+Display:ital,wght@0,400;0,700;1,400&display=swap'
          }
        ]
      }
  }]
}