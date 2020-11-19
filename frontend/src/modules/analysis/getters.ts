import { Getters } from "vuex-smart-module";
import { AnalysisState } from ".";

export class AnalysisGetters extends Getters<AnalysisState> {
  get regionsEnabled() {
    return this.state.regionsEnabled;
  }

  get selectedRegionStats() {
    return this.state.selectedRegionStats;
  }
}
