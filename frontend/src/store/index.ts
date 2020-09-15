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
import { responsiveModule } from "@/modules/responsive";
import { RootActions } from "@/store/actions";
import { selectionModule } from "@/modules/selection";
import { centroidsModule } from "@/modules/centroids";
import { presetModule } from "@/modules/presets";
import { gateModule } from "@/modules/gates";
import { groupModule } from "@/modules/group";
import { memberModule } from "@/modules/member";
import { resultModule } from "@/modules/results";

Vue.use(Vuex);

const debug = false; // process.env.NODE_ENV !== "production";

const rootModule = new Module({
  modules: {
    group: groupModule,
    member: memberModule,
    analysis: analysisModule,
    dataset: datasetModule,
    results: resultModule,
    experiment: experimentModule,
    main: mainModule,
    responsive: responsiveModule,
    settings: settingsModule,
    user: userModule,
    selection: selectionModule,
    centroids: centroidsModule,
    presets: presetModule,
    gates: gateModule,
  },
  actions: RootActions,
});

const vuexStorage = new VuexPersistence<typeof rootModule>({
  strictMode: debug,
  storage: localforage as any,
  asyncStorage: true,
  modules: ["settings"],
});

export const store = createStore(rootModule, {
  strict: debug,
  plugins: debug ? [vuexStorage.plugin, createLogger()] : [vuexStorage.plugin],
  mutations: {
    RESTORE_MUTATION: vuexStorage.RESTORE_MUTATION, // this mutation **MUST** be named "RESTORE_MUTATION"
  },
});
export default store;
