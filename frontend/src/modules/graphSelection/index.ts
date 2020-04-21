import { Module } from "vuex-smart-module";
import { GraphSelectionActions } from "./actions";
import { GraphSelectionGetters } from "./getters";
import { GraphSelectionMutations } from "./mutations";

export class GraphSelectionState {
  tool = "lasso"; // what selection tool mode (lasso, brush, ...)
  selection: any = { mode: "all" }; // current selection, which is tool specific
}

export const graphSelectionModule = new Module({
  namespaced: true,

  state: GraphSelectionState,
  getters: GraphSelectionGetters,
  mutations: GraphSelectionMutations,
  actions: GraphSelectionActions,
});
