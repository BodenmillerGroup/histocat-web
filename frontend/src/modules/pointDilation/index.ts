import { Module } from "vuex-smart-module";
import { PointDilationActions } from "./actions";
import { PointDilationGetters } from "./getters";
import { PointDilationMutations } from "./mutations";


export class PointDilationState {
  metadataField: string = "";
  categoryField: string = "";
}

export const pointDilationModule = new Module({
  namespaced: true,

  state: PointDilationState,
  getters: PointDilationGetters,
  mutations: PointDilationMutations,
  actions: PointDilationActions,
});
