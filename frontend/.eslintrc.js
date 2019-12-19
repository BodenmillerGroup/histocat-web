module.exports = {
  env: {
    "browser": true,
    "commonjs": true,
    "jest": true,
    "serviceworker": true
  },
  parser: "vue-eslint-parser",
  parserOptions: {
    parser: "@typescript-eslint/parser",
    ecmaVersion: 2018, // Allows for the parsing of modern ECMAScript features
    sourceType: "module", // Allows for the use of imports
    extraFileExtensions: [
      ".vue"
    ]
  }
};
