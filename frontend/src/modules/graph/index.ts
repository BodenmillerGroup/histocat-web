import { Module } from "vuex-smart-module";
import { GraphActions } from "./actions";
import { GraphGetters } from "./getters";
import { GraphMutations } from "./mutations";

export class GraphState {

}

export const graphModule = new Module({
  namespaced: true,

  state: GraphState,
  getters: GraphGetters,
  mutations: GraphMutations,
  actions: GraphActions,
});
