import { Module } from "vuex-smart-module";
import { ResultActions } from "./actions";
import { ResultGetters } from "./getters";
import { IResult } from "./models";
import { ResultMutations } from "./mutations";
import { schema } from "normalizr";

export const resultSchema = new schema.Entity("results");
export const resultListSchema = [resultSchema];

export class ResultState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IResult } = {};
  activeResultId: number | null = null;
}

export const resultModule = new Module({
  namespaced: true,

  state: ResultState,
  getters: ResultGetters,
  mutations: ResultMutations,
  actions: ResultActions,
});
