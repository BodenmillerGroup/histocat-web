import { Module } from "vuex-smart-module";
import { ExperimentActions } from "./actions";
import { ExperimentGetters } from "./getters";
import { IExperiment } from "./models";
import { ExperimentMutations } from "./mutations";

export class ExperimentState {
  experiments: IExperiment[] = [];
  tags: string[] = [];

  activeExperimentId?: number = undefined;
  activeAcquisitionId?: number = undefined;
  activeWorkspaceNode?: object = undefined;
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
