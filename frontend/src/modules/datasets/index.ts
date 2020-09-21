import { Module } from "vuex-smart-module";
import { DatasetsActions } from "./actions";
import { DatasetsGetters } from "./getters";
import { IDataset } from "./models";
import { DatasetsMutations } from "./mutations";
import { schema } from "normalizr";

export const datasetSchema = new schema.Entity("datasets");
export const datasetListSchema = [datasetSchema];

export class DatasetsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IDataset } = {};
  activeDatasetId: number | null = null;
}

export const datasetsModule = new Module({
  namespaced: true,

  state: DatasetsState,
  getters: DatasetsGetters,
  mutations: DatasetsMutations,
  actions: DatasetsActions,
});
