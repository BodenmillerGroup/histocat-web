import { Module } from "vuex-smart-module";
import { PipelinesActions } from "./actions";
import { PipelinesGetters } from "./getters";
import { IPipeline } from "./models";
import { PipelinesMutations } from "./mutations";
import { schema } from "normalizr";

export const pipelineSchema = new schema.Entity("pipelines");
export const pipelineListSchema = [pipelineSchema];

export class PipelinesState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IPipeline } = {};
  activePipelineId: number | null = null;

  selectedAcquisitionIds: number[] = [];
  steps: any[] = [];
}

export const pipelinesModule = new Module({
  namespaced: true,

  state: PipelinesState,
  getters: PipelinesGetters,
  mutations: PipelinesMutations,
  actions: PipelinesActions,
});
