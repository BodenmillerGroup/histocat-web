import { Module } from "vuex-smart-module";
import { DatasetActions } from "./actions";
import { DatasetGetters } from "./getters";
import { IDataset } from "./models";
import { DatasetMutations } from "./mutations";

export class DatasetState {
  datasets: IDataset[] = [];
  activeDataset?: IDataset = undefined;
}

export const datasetModule = new Module({
  namespaced: false,

  state: DatasetState,
  getters: DatasetGetters,
  mutations: DatasetMutations,
  actions: DatasetActions
});
