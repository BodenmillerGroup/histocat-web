import { Getters } from "vuex-smart-module";
import { AnalysisState } from ".";

export class AnalysisGetters extends Getters<AnalysisState> {
  get scatterPlotData() {
    return this.state.scatterPlotData;
  }

  get boxPlotData() {
    return this.state.boxPlotData;
  }

  get pcaData() {
    return this.state.pcaData;
  }

  get tsneData() {
    return this.state.tsneData;
  }

  get umapData() {
    return this.state.umapData;
  }

  get phenographData() {
    return this.state.phenographData;
  }

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
