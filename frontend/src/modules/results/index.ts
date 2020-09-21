import { Module } from "vuex-smart-module";
import { ResultsActions } from "./actions";
import { ResultsGetters } from "./getters";
import { IResult } from "./models";
import { ResultsMutations } from "./mutations";
import { schema } from "normalizr";

export const resultSchema = new schema.Entity("results");
export const resultListSchema = [resultSchema];

export class ResultsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IResult } = {};
  activeResultId: number | null = null;
}

export const resultsModule = new Module({
  namespaced: true,

  state: ResultsState,
  getters: ResultsGetters,
  mutations: ResultsMutations,
  actions: ResultsActions,
});
