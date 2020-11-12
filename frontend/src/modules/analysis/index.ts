import Feature from "ol/Feature";
import { Module } from "vuex-smart-module";
import { AnalysisActions } from "./actions";
import { AnalysisGetters } from "./getters";
import {
  IPhenoGraphData,
  IPlotSeries,
  IRegionChannelData,
  IScatterPlotData,
} from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisState {
  scatterPlotData: IScatterPlotData | null = null;
  boxPlotData: IPlotSeries[] = [];
  phenographData: IPhenoGraphData | null = null;

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
