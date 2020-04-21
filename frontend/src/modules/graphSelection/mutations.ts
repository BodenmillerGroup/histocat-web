import { Mutations } from "vuex-smart-module";
import { GraphSelectionState } from ".";

export class GraphSelectionMutations extends Mutations<GraphSelectionState> {
  setTool(value: string) {
    this.state.tool = value;
  }

  setSelection(value: any) {
    this.state.selection = value;
  }
}
