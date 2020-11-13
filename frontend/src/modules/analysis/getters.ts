import { Getters } from "vuex-smart-module";
import { AnalysisState } from ".";

export class AnalysisGetters extends Getters<AnalysisState> {
  get regionsEnabled() {
    return this.state.regionsEnabled;
  }

  get selectedRegion() {
    return this.state.selectedRegion;
  }

  get selectedRegionStats() {
    return this.state.selectedRegionStats;
  }
}
