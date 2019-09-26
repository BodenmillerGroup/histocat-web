import Feature from "ol/Feature";
import { Module } from "vuex-smart-module";
import { AnalysisActions } from "./actions";
import { AnalysisGetters } from "./getters";
import { IPCAData, IPlotSeries, IRegionChannelData, IScatterPlotData, ITSNEData, IUMAPData } from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisState {
  segmentationImage: string | ArrayBuffer | null = null;
  segmentationContours: number[][] = [];
  scatterPlotData: IScatterPlotData | null = null;
  boxPlotData: IPlotSeries[] = [];
  pcaData: IPCAData | null = null;
  tsneData: ITSNEData | null = null;
  umapData: IUMAPData | null = null;

  regionsEnabled = false;
  selectedRegion: Feature | null = null;
  selectedRegionStats: IRegionChannelData[] = [];
}

export const analysisModule = new Module({
  namespaced: false,

  state: AnalysisState,
  getters: AnalysisGetters,
  mutations: AnalysisMutations,
  actions: AnalysisActions
});
