import { Module } from "vuex-smart-module";
import { DatasetActions } from "./actions";
import { DatasetGetters } from "./getters";
import { IDataset } from "./models";
import { DatasetMutations } from "./mutations";
import { schema } from "normalizr";

export const datasetSchema = new schema.Entity("datasets");
export const datasetListSchema = [datasetSchema];

export class DatasetState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IDataset } = {};
  activeDatasetId: number | null = null;
}

export const datasetModule = new Module({
  namespaced: true,

  state: DatasetState,
  getters: DatasetGetters,
  mutations: DatasetMutations,
  actions: DatasetActions,
});
