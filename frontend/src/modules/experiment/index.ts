import { Module } from "vuex-smart-module";
import { ExperimentActions } from "./actions";
import { ExperimentGetters } from "./getters";
import { IExperiment, IExperimentData } from "./models";
import { ExperimentMutations } from "./mutations";
import { schema } from "normalizr";

export const experimentSchema = new schema.Entity("experiments");
export const experimentListSchema = [experimentSchema];

export class ExperimentState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IExperiment } = {};
  tags: string[] = [];
  experimentData: IExperimentData | null = null;

  activeExperimentId: number | null = null;
  activeAcquisitionId: number | null = null;
  activeWorkspaceNode: object | null = null;
  selectedAcquisitionIds: number[] = [];
  selectedMetals: string[] = [];
  channelStackImage: string | ArrayBuffer | null = null;

  colorizeMaskInProgress = false;
}

export const experimentModule = new Module({
  namespaced: true,

  state: ExperimentState,
  getters: ExperimentGetters,
  mutations: ExperimentMutations,
  actions: ExperimentActions,
});
