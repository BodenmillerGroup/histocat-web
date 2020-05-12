import Feature from "ol/Feature";
import { Mutations } from "vuex-smart-module";
import { AnalysisState } from ".";
import {
  ICentroidsData,
  IPCAData,
  IPhenoGraphData,
  IPlotSeries,
  IRegionChannelData,
  IScatterPlotData,
  ITSNEData,
  IUMAPData,
} from "./models";

export class AnalysisMutations extends Mutations<AnalysisState> {
  setCentroidsData(data: ICentroidsData | null) {
    this.state.centroidsData = data;
  }

  setSegmentationImage(base64Image: string | ArrayBuffer | null) {
    this.state.segmentationImage = base64Image;
  }

  setSegmentationContours(contours: number[][]) {
    this.state.segmentationContours = contours;
  }

  setScatterPlotData(data: IScatterPlotData | null) {
    this.state.scatterPlotData = data;
  }

  setBoxPlotData(data: IPlotSeries[]) {
    this.state.boxPlotData = data;
  }

  setPCAData(data: IPCAData | null) {
    this.state.pcaData = data;
  }

  setTSNEData(data: ITSNEData | null) {
    this.state.tsneData = data;
  }

  setUMAPData(data: IUMAPData | null) {
    this.state.umapData = data;
  }

  setPhenoGraphData(data: IPhenoGraphData | null) {
    this.state.phenographData = data;
  }

  setRegionsEnabled(state: boolean) {
    this.state.regionsEnabled = state;
  }

  setSelectedRegion(region: Feature | null) {
    this.state.selectedRegion = region;
  }

  setSelectedRegionStats(stats: IRegionChannelData[]) {
    this.state.selectedRegionStats = stats;
  }

  reset() {
    this.state.segmentationImage = null;
    this.state.segmentationContours = [];
    this.state.scatterPlotData = null;
    this.state.boxPlotData = [];
    this.state.pcaData = null;
    this.state.tsneData = null;
    this.state.umapData = null;
    this.state.phenographData = null;
    this.state.regionsEnabled = false;
    this.state.selectedRegion = null;
    this.state.selectedRegionStats = [];
  }
}
