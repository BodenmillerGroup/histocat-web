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
  selectedExperimentId?: number;
  selectedAcquisitionId?: number;
  selectedMetals: string[];
}

const defaultState: ExperimentsState = {
  experiments: [],
  datasets: [],
  tags: [],
  selectedExperimentId: undefined,
  selectedAcquisitionId: undefined,
  selectedMetals: [],
};

export const experimentModule: Module<ExperimentsState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
