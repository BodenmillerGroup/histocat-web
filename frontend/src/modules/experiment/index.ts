import { Module } from 'vuex-smart-module';
import { ExperimentActions } from './actions';
import { ExperimentGetters } from './getters';
import { IDataset, IExperiment, IShare } from './models';
import { ExperimentMutations } from './mutations';

export class ExperimentState {
  experiments: IExperiment[] = [];
  datasets: IDataset[] = [];
  tags: string[] = [];
  shares: IShare[] = [];
  activeExperimentId?: number = undefined;
  activeSlideId?: number = undefined;
  activePanoramaId?: number = undefined;
  activeAcquisitionId?: number = undefined;
  activeWorkspaceNode?: object = undefined;
  selectedAcquisitionIds: number[] = [];
  selectedMetals: string[] = [];
  channelStackImage: string | ArrayBuffer | null = null;
}

export const experimentModule = new Module({
  namespaced: false,

  state: ExperimentState,
  getters: ExperimentGetters,
  mutations: ExperimentMutations,
  actions: ExperimentActions,
});
