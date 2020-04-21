import { Module } from "vuex-smart-module";
import { CrossfilterActions } from "./actions";
import { CrossfilterGetters } from "./getters";
import { CrossfilterMutations } from "./mutations";

export class CrossfilterState {

}

export const crossfilterModule = new Module({
  namespaced: true,

  state: CrossfilterState,
  getters: CrossfilterGetters,
  mutations: CrossfilterMutations,
  actions: CrossfilterActions,
});
