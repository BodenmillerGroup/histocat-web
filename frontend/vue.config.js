module.exports = {
  lintOnSave: false,
  runtimeCompiler: false,
  transpileDependencies: [
    "vuetify"
  ],

  pwa: {
    name: "HistoCAT"
  },

  configureWebpack: config => {
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
        source: false
      }
    }
  },

  devServer: {
    port: 9999
  }
};
