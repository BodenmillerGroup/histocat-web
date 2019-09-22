import { Module } from "vuex-smart-module";
import { ExperimentActions } from "./actions";
import { ExperimentGetters } from "./getters";
import { IExperiment, IShare } from "./models";
import { ExperimentMutations } from "./mutations";
import Feature from 'ol/Feature';

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

  features: Map<string, Feature> = new Map<string, Feature>();
}

export const experimentModule = new Module({
  namespaced: false,

  state: ExperimentState,
  getters: ExperimentGetters,
  mutations: ExperimentMutations,
  actions: ExperimentActions
});
