import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { SelectionState } from ".";
import { SelectionGetters } from "./getters";
import { SelectionMutations } from "./mutations";
import { SelectedCell } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SELECTED_CELLS } from "./events";

export class SelectionActions extends Actions<SelectionState, SelectionGetters, SelectionMutations, SelectionActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}

  setSelectedCells(payload: SelectedCell[], isGlobal = true) {
    BroadcastManager.publish(SET_SELECTED_CELLS, payload, isGlobal);
  }
}
