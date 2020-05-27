import { Getters } from "vuex-smart-module";
import { SelectionState } from ".";

export class SelectionGetters extends Getters<SelectionState> {
  get cellIds() {
    return this.state.cellIds;
  }
}
