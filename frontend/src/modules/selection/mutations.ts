import { Mutations } from "vuex-smart-module";
import { SelectionState } from ".";
import { SelectedCell } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SELECTED_CELLS } from "./events";

export class SelectionMutations extends Mutations<SelectionState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_SELECTED_CELLS, (payload) => this.setSelectedCells(payload));
  }

  setSelectedCells(payload: SelectedCell[]) {
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
