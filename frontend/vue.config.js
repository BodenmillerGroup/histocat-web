const manifestJSON = require("./public/manifest.json");

module.exports = {
  lintOnSave: false,
  runtimeCompiler: false,
  transpileDependencies: ["vuetify"],

  pwa: {
    themeColor: manifestJSON.theme_color,
    name: manifestJSON.short_name,
    msTileColor: manifestJSON.background_color,
    appleMobileWebAppCapable: "yes",
    appleMobileWebAppStatusBarStyle: "black",
    workboxPluginMode: "InjectManifest",
    workboxOptions: {
      swSrc: "public/service-worker.js",
    },
  },

  configureWebpack: (config) => {
    if (process.env.NODE_ENV === "production") {
      ("nosources-source-map");
    } else {
      ("eval-source-map");
    }
  },

  pluginOptions: {
    webpackBundleAnalyzer: {
      analyzerMode: "disabled",
      generateStatsFile: false,
      openAnalyzer: false,
      // Excludes module sources from stats file so there won't be any sensitive data
      statsOptions: {
        source: false,
      },
    },
  },

  devServer: {
    port: 9999,
  },
};
