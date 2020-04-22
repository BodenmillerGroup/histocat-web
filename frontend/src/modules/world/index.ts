import { Module } from "vuex-smart-module";
import { WorldActions } from "./actions";
import { WorldGetters } from "./getters";
import { WorldMutations } from "./mutations";
import {IWorld} from "@/modules/world/models";

export class WorldState {
  world: IWorld | null = null;
}

export const worldModule = new Module({
  namespaced: true,

  state: WorldState,
  getters: WorldGetters,
  mutations: WorldMutations,
  actions: WorldActions,
});
