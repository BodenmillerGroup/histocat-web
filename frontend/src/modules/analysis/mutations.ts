import { Mutations } from 'vuex-smart-module';
import { AnalysisState } from '.';
import { IScatterPlotData } from './models';


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

  resetAnalysis() {
    this.setSegmentationImage(null);
    this.setScatterPlotData(null);
  }
}
