import { Module } from "vuex-smart-module";
import { WorldActions } from "./actions";
import { WorldGetters } from "./getters";
import { WorldMutations } from "./mutations";

export class WorldState {

}

export const worldModule = new Module({
  namespaced: true,

  state: WorldState,
  getters: WorldGetters,
  mutations: WorldMutations,
  actions: WorldActions,
});
