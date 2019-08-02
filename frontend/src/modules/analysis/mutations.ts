import { Mutations } from 'vuex-smart-module';
import { AnalysisState } from '.';


export class AnalysisMutations extends Mutations<AnalysisState> {
  setAnalysisImage(base64Image: string | ArrayBuffer | null) {
    this.state.analysisImage = base64Image;
  }

  resetAnalysis() {
    this.setAnalysisImage(null);
  }
}
