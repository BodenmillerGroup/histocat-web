module.exports = {
  preset: "@vue/cli-plugin-unit-jest/presets/typescript-and-babel",
  moduleNameMapper: {
    "^ky$": require.resolve("ky").replace("index.js", "umd.js"),
  },
};
