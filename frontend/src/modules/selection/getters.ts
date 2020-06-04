import { Getters } from "vuex-smart-module";
import { SelectionState } from ".";

export class SelectionGetters extends Getters<SelectionState> {
  get selectedCells() {
    return this.state.selectedCells;
  }
}
