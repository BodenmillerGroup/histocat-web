import { Module } from "vuex-smart-module";
import { AnalysisActions } from "./actions";
import { AnalysisGetters } from "./getters";
import { IRegionChannelData } from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisState {
  regionsEnabled = false;
  selectedRegionStats: IRegionChannelData[] = [];
}

export const analysisModule = new Module({
  namespaced: true,

  state: AnalysisState,
  getters: AnalysisGetters,
  mutations: AnalysisMutations,
  actions: AnalysisActions,
});
