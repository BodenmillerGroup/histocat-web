import { Mutations } from 'vuex-smart-module';
import { AnalysisState } from '.';
import { IPCAData, IPlotSeries, IScatterPlotData } from './models';


export class AnalysisMutations extends Mutations<AnalysisState> {
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

  reset() {
    this.setSegmentationImage(null);
    this.setSegmentationContours([]);
    this.setScatterPlotData(null);
    this.setBoxPlotData([]);
    this.setPCAData(null);
  }
}
