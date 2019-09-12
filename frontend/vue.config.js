module.exports = {
  assetsDir: 'assets',
  runtimeCompiler: false,

  pwa: {
    name: 'HistoCAT',
  },

  configureWebpack: config => {
    if (process.env.NODE_ENV === 'production') {
      devtool: 'nosources-source-map';
    } else {
      devtool: 'eval-source-map';
    }
  },

  pluginOptions: {
    webpackBundleAnalyzer: {
      openAnalyzer: false,
    },
  },
};
