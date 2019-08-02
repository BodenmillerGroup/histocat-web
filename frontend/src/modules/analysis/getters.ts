import { Getters } from 'vuex-smart-module';
import { AnalysisState } from '.';

export class AnalysisGetters extends Getters<AnalysisState> {
  get analysisImage() {
    return this.state.analysisImage;
  }
}
