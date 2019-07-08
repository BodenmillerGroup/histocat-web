import { IDataset, IExperiment } from '@/modules/experiment/models';
import { RootState } from '@/store';
import { Module } from 'vuex';
import { actions } from './actions';
import { getters } from './getters';
import { mutations } from './mutations';

export interface ExperimentsState {
  experiments: IExperiment[];
  datasets: IDataset[];
  tags: string[];
  activeExperimentId?: number;
  activeSlideId?: number;
  activePanoramaId?: number;
  activeAcquisitionId?: number;
  selectedAcquisitionIds: number[];
  selectedMetals: string[];
}

const defaultState: ExperimentsState = {
  experiments: [],
  datasets: [],
  tags: [],
  activeExperimentId: undefined,
  activeSlideId: undefined,
  activePanoramaId: undefined,
  activeAcquisitionId: undefined,
  selectedAcquisitionIds: [],
  selectedMetals: [],
};

export const experimentModule: Module<ExperimentsState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
