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
      openAnalyzer: false
    }
  },

  devServer: {
    port: 9999
  }
};
