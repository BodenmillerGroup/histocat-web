import Feature from "ol/Feature";
import { Module } from "vuex-smart-module";
import { AnalysisActions } from "./actions";
import { AnalysisGetters } from "./getters";
import { IRegionChannelData } from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisState {
  regionsEnabled = false;
  selectedRegion: Feature | null = null;
  selectedRegionStats: IRegionChannelData[] = [];
}

export const analysisModule = new Module({
  namespaced: true,

  state: AnalysisState,
  getters: AnalysisGetters,
  mutations: AnalysisMutations,
  actions: AnalysisActions,
});
