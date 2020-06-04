import { Module } from "vuex-smart-module";
import { SelectionActions } from "./actions";
import { SelectionGetters } from "./getters";
import { SelectionMutations } from "./mutations";
import { SelectedCell } from "./models";

export class SelectionState {
  // Mapped by acquisition id
  selectedCells: Map<number, SelectedCell[]> | null = null;
}

export const selectionModule = new Module({
  namespaced: true,

  state: SelectionState,
  getters: SelectionGetters,
  mutations: SelectionMutations,
  actions: SelectionActions,
});
