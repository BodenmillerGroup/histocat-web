module.exports = {
  preset: "@vue/cli-plugin-unit-jest/presets/typescript-and-babel",
  transformIgnorePatterns: ["/node_modules/(?!ky|lodash-es)"],
  setupFiles: ["<rootDir>/tests/jest.stubs.js"],
  snapshotSerializers: ["jest-serializer-vue"],
};
