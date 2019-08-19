import { Mutations } from 'vuex-smart-module';
import { AnalysisState } from '.';


export class AnalysisMutations extends Mutations<AnalysisState> {
  setSegmentationImage(base64Image: string | ArrayBuffer | null) {
    this.state.segmentationImage = base64Image;
  }

  setSegmentationContours(contours: number[][]) {
    this.state.segmentationContours = contours;
  }

  resetAnalysis() {
    this.setSegmentationImage(null);
  }
}
