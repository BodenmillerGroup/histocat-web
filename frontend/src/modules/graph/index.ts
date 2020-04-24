import { Module } from "vuex-smart-module";
import { GraphActions } from "./actions";
import { GraphGetters } from "./getters";
import { GraphMutations } from "./mutations";
import Dataframe from "@/cellxgene/util/dataframe/dataframe";

export class GraphState {
  schema: any = null;
  obsAnnotations: Dataframe | null = null;
  varAnnotations: Dataframe | null = null;
}

export const graphModule = new Module({
  namespaced: true,

  state: GraphState,
  getters: GraphGetters,
  mutations: GraphMutations,
  actions: GraphActions,
});
