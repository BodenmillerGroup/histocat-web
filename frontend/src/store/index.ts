import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { userModule } from "@/modules/user";
import Vue from "vue";
import Vuex from "vuex";
import VuexPersistence from "vuex-persist";
import { createStore, Module } from "vuex-smart-module";
import createLogger from "vuex/dist/logger";
import { responsiveModule } from "@/modules/responsive";
import { RootActions } from "@/store/actions";
import { selectionModule } from "@/modules/selection";
import { centroidsModule } from "@/modules/centroids";
import { presetsModule } from "@/modules/presets";
import { gatesModule } from "@/modules/gates";
import { groupModule } from "@/modules/group";
import { memberModule } from "@/modules/member";
import { resultsModule } from "@/modules/results";

Vue.use(Vuex);

const debug = false; // process.env.NODE_ENV !== "production";

const rootModule = new Module({
  modules: {
    group: groupModule,
    member: memberModule,
    analysis: analysisModule,
    dataset: datasetsModule,
    results: resultsModule,
    projects: projectsModule,
    main: mainModule,
    responsive: responsiveModule,
    settings: settingsModule,
    user: userModule,
    selection: selectionModule,
    centroids: centroidsModule,
    presets: presetsModule,
    gates: gatesModule,
  },
  actions: RootActions,
});

const vuexStorage = new VuexPersistence<typeof rootModule>({
  strictMode: debug,
  storage: window.localStorage,
  asyncStorage: false,
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
