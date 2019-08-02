import { Module } from 'vuex-smart-module';
import { AnalysisActions } from './actions';
import { AnalysisGetters } from './getters';
import { AnalysisMutations } from './mutations';

export class AnalysisState {
  analysisImage: string | ArrayBuffer | null = null;
}

export const analysisModule = new Module({
  namespaced: false,

  state: AnalysisState,
  getters: AnalysisGetters,
  mutations: AnalysisMutations,
  actions: AnalysisActions,
});
