import { Module } from "vuex-smart-module";
import { CentroidsActions } from "./actions";
import { CentroidsGetters } from "./getters";
import { CentroidsMutations } from "./mutations";
import { CellPoint } from "@/data/CellPoint";

export class CentroidsState {
  // Map centroids by acquisitionId
  centroids: Map<number, CellPoint[]> | null = null;
}

export const centroidsModule = new Module({
  namespaced: true,

  state: CentroidsState,
  getters: CentroidsGetters,
  mutations: CentroidsMutations,
  actions: CentroidsActions,
});
