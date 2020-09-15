import Feature from "ol/Feature";
import { Module } from "vuex-smart-module";
import { AnalysisActions } from "./actions";
import { AnalysisGetters } from "./getters";
import {
  IPCAData,
  IPhenoGraphData,
  IPlotSeries,
  IRegionChannelData,
  IScatterPlotData,
  ITSNEData,
  IUMAPData,
} from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisState {
  scatterPlotData: IScatterPlotData | null = null;
  boxPlotData: IPlotSeries[] = [];
  pcaData: IPCAData | null = null;
  tsneData: ITSNEData | null = null;
  umapData: IUMAPData | null = null;
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
