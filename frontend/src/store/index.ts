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
import { graphModule } from "@/modules/graph";
import { responsiveModule } from "@/modules/responsive";
import { centroidLabelsModule } from "@/modules/centroidLabels";
import { colorsModule } from "@/modules/colors";
import { pointDilationModule } from "@/modules/pointDilation";
import { worldModule } from "@/modules/world";
import { universeModule } from "@/modules/universe";
import { graphSelectionModule } from "@/modules/graphSelection";
import { layoutChoiceModule } from "@/modules/layoutChoice";
import { controlsModule } from "@/modules/controls";
import { crossfilterModule } from "@/modules/crossfilter";

Vue.use(Vuex);

const debug = false; // process.env.NODE_ENV !== "production";

const rootModule = new Module({
  modules: {
    analysis: analysisModule,
    centroidLabels: centroidLabelsModule,
    colors: colorsModule,
    controls: controlsModule,
    crossfilter: crossfilterModule,
    dataset: datasetModule,
    experiment: experimentModule,
    graph: graphModule,
    graphSelection: graphSelectionModule,
    layoutChoice: layoutChoiceModule,
    main: mainModule,
    pointDilation: pointDilationModule,
    responsive: responsiveModule,
    settings: settingsModule,
    universe: universeModule,
    user: userModule,
    world: worldModule,
  },
});

const vuexStorage = new VuexPersistence<typeof rootModule>({
  strictMode: debug,
  storage: localforage,
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
