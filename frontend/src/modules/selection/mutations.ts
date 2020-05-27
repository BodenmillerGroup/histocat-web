import { Mutations } from "vuex-smart-module";
import { SelectionState } from ".";

export class SelectionMutations extends Mutations<SelectionState> {
  setCellIds(payload: Map<number, number[]>) {
    console.log(payload)
    this.state.cellIds = payload;
  }

  reset() {
    // acquire initial state
    const s = new SelectionState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
