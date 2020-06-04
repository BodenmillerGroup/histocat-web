import { Mutations } from "vuex-smart-module";
import { SelectionState } from ".";
import { SelectedCell } from "./models";

export class SelectionMutations extends Mutations<SelectionState> {
  setSelectedCells(payload: Map<number, SelectedCell[]> | null) {
    this.state.selectedCells = payload;
  }

  reset() {
    // acquire initial state
    const s = new SelectionState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
