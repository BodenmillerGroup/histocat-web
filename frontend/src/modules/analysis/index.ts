import { Module } from 'vuex-smart-module';
import { AnalysisActions } from './actions';
import { AnalysisGetters } from './getters';
import { IPlotSeries, IScatterPlotData } from './models';
import { AnalysisMutations } from './mutations';

export class AnalysisState {
  segmentationImage: string | ArrayBuffer | null = null;
  segmentationContours: number[][] = [];
  scatterPlotData: IScatterPlotData | null = null;
  boxPlotData: IPlotSeries[] = [];
}

export const analysisModule = new Module({
  namespaced: false,

  state: AnalysisState,
  getters: AnalysisGetters,
  mutations: AnalysisMutations,
  actions: AnalysisActions,
});
