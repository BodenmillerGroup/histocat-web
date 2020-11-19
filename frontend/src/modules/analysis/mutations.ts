import { Mutations } from "vuex-smart-module";
import { AnalysisState } from ".";
import { IRegionChannelData } from "./models";

export class AnalysisMutations extends Mutations<AnalysisState> {
  setRegionsEnabled(state: boolean) {
    this.state.regionsEnabled = state;
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
