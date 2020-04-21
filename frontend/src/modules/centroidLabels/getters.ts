import { Getters } from "vuex-smart-module";
import { CentroidLabelsState } from ".";

export class CentroidLabelsGetters extends Getters<CentroidLabelsState> {
  get labels() {
    return this.state.labels;
  }

  get showLabels() {
    return this.state.showLabels;
  }
}
