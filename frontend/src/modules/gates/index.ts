import { Module } from "vuex-smart-module";
import { GateActions } from "./actions";
import { GateGetters } from "./getters";
import { IGate } from "./models";
import { GateMutations } from "./mutations";
import { schema } from "normalizr";

export const gateSchema = new schema.Entity("gates");
export const gateListSchema = [gateSchema];

export class GateState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IGate } = {};
  activeGateId: number | null = null;
}

export const gateModule = new Module({
  namespaced: true,

  state: GateState,
  getters: GateGetters,
  mutations: GateMutations,
  actions: GateActions,
});
