import { Module } from "vuex-smart-module";
import { CentroidLabelsActions } from "./actions";
import { CentroidLabelsGetters } from "./getters";
import { CentroidLabelsMutations } from "./mutations";

export class CentroidLabelsState {
  labels: Map<any, any> = new Map<any, any>();
  showLabels = false;
}

export const centroidLabelsModule = new Module({
  namespaced: true,

  state: CentroidLabelsState,
  getters: CentroidLabelsGetters,
  mutations: CentroidLabelsMutations,
  actions: CentroidLabelsActions,
});
