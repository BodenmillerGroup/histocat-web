import { Module } from "vuex-smart-module";
import { UniverseActions } from "./actions";
import { UniverseGetters } from "./getters";
import { UniverseMutations } from "./mutations";
import { IUniverse } from "@/modules/universe/models";

export class UniverseState {
  universe: IUniverse | null = null;
}

export const universeModule = new Module({
  namespaced: true,

  state: UniverseState,
  getters: UniverseGetters,
  mutations: UniverseMutations,
  actions: UniverseActions,
});
