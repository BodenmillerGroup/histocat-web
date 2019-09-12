import { Module } from "vuex-smart-module";
import { ExperimentActions } from "./actions";
import { ExperimentGetters } from "./getters";
import { IExperiment, IShare } from "./models";
import { ExperimentMutations } from "./mutations";

export class ExperimentState {
  experiments: IExperiment[] = [];
  tags: string[] = [];
  shares: IShare[] = [];

  activeExperimentId?: number = undefined;
  activeAcquisitionId?: number = undefined;
  activeWorkspaceNode?: object = undefined;
  selectedAcquisitionIds: number[] = [];
  selectedMetals: string[] = [];
  channelStackImage: string | ArrayBuffer | null = null;

  colorizeMaskInProgress = false;
}

export const experimentModule = new Module({
  namespaced: false,

  state: ExperimentState,
  getters: ExperimentGetters,
  mutations: ExperimentMutations,
  actions: ExperimentActions
});
