import { Module } from "vuex-smart-module";
import { GatesActions } from "./actions";
import { GatesGetters } from "./getters";
import { IGate } from "./models";
import { GatesMutations } from "./mutations";
import { schema } from "normalizr";

export const gateSchema = new schema.Entity("gates");
export const gateListSchema = [gateSchema];

export class GatesState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IGate } = {};
  activeGateId: number | null = null;
}

export const gatesModule = new Module({
  namespaced: true,

  state: GatesState,
  getters: GatesGetters,
  mutations: GatesMutations,
  actions: GatesActions,
});
