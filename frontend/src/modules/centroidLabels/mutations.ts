import { Mutations } from "vuex-smart-module";
import { CentroidLabelsState } from ".";

export class CentroidLabelsMutations extends Mutations<CentroidLabelsState> {
  setLabels(labels: Map<any, any>) {
    this.state.labels = labels;
  }

  setShowLabels(state: boolean) {
    this.state.showLabels = state;
  }

  resetLabels() {
    this.state.labels = new Map<any, any>();
  }

  reset() {
    this.state.labels = new Map<any, any>();
    this.state.showLabels = false;
  }
}
