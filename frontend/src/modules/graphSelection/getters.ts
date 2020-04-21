import { Getters } from "vuex-smart-module";
import { GraphSelectionState } from ".";

export class GraphSelectionGetters extends Getters<GraphSelectionState> {
  get tool() {
    return this.state.tool;
  }

  get selection() {
    return this.state.selection;
  }
}
