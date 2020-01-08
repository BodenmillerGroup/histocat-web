import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { userModule } from "@/modules/user";
import localforage from "localforage";
import Vue from "vue";
import Vuex from "vuex";
import VuexPersistence from "vuex-persist";
import { createStore, Module } from "vuex-smart-module";
import createLogger from "vuex/dist/logger";

Vue.use(Vuex);

const debug = false; // process.env.NODE_ENV !== "production";

const rootModule = new Module({
  modules: {
    main: mainModule,
    user: userModule,
    experiment: experimentModule,
    dataset: datasetModule,
    analysis: analysisModule,
    settings: settingsModule
  }
});

const vuexStorage = new VuexPersistence<typeof rootModule>({
  strictMode: debug,
  storage: localforage,
  asyncStorage: true,
  modules: ["settings"]
});

export const store = createStore(rootModule, {
  strict: debug,
  plugins: debug ? [vuexStorage.plugin, createLogger()] : [vuexStorage.plugin],
  mutations: {
    RESTORE_MUTATION: vuexStorage.RESTORE_MUTATION // this mutation **MUST** be named "RESTORE_MUTATION"
  }
});
export default store;
