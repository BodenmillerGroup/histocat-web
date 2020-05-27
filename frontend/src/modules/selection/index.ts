import { Module } from "vuex-smart-module";
import { SelectionActions } from "./actions";
import { SelectionGetters } from "./getters";
import { SelectionMutations } from "./mutations";

export class SelectionState {
  cellIds: Map<number, number[]> | null = null;
}

export const selectionModule = new Module({
  namespaced: true,

  state: SelectionState,
  getters: SelectionGetters,
  mutations: SelectionMutations,
  actions: SelectionActions,
});
