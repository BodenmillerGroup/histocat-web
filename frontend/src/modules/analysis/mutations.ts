import Feature from "ol/Feature";
import { Mutations } from "vuex-smart-module";
import { AnalysisState } from ".";
import {
  IPCAData,
  IPhenoGraphData,
  IPlotSeries,
  IRegionChannelData,
  IScatterPlotData,
  ITSNEData,
  IUMAPData,
} from "./models";

export class AnalysisMutations extends Mutations<AnalysisState> {
  setSegmentationImage(base64Image: string | ArrayBuffer | null) {
    this.state.segmentationImage = base64Image;
  }

  setSegmentationContours(contours: number[][]) {
    this.state.segmentationContours = contours;
  }

  setScatterPlotData(data: IScatterPlotData | null) {
    this.state.scatterPlotData = Object.freeze(data);
  }

  setBoxPlotData(data: IPlotSeries[]) {
    this.state.boxPlotData = data;
  }

  setPCAData(data: IPCAData | null) {
    this.state.pcaData = Object.freeze(data);
  }

  setTSNEData(data: ITSNEData | null) {
    this.state.tsneData = Object.freeze(data);
  }

  setUMAPData(data: IUMAPData | null) {
    this.state.umapData = Object.freeze(data);
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
    // acquire initial state
    const s = new AnalysisState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
