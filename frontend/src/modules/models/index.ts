import { Module } from "vuex-smart-module";
import { ModelsActions } from "./actions";
import { ModelsGetters } from "./getters";
import { IModel } from "./models";
import { ModelsMutations } from "./mutations";
import { schema } from "normalizr";

export const modelSchema = new schema.Entity("models");
export const modelListSchema = [modelSchema];

export class ModelsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IModel } = {};
  activeModelId: number | null = null;
}

export const modelsModule = new Module({
  namespaced: true,

  state: ModelsState,
  getters: ModelsGetters,
  mutations: ModelsMutations,
  actions: ModelsActions,
});
