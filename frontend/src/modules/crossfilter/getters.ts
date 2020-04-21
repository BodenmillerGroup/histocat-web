import { Getters } from "vuex-smart-module";
import { CrossfilterState } from ".";

export class CrossfilterGetters extends Getters<CrossfilterState> {
  get crossfilter() {
    return this.state.crossfilter;
  }
}
