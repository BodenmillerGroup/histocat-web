import { mutations } from './mutations';
import { getters } from './getters';
import { actions } from './actions';
import { ExperimentsState } from './state';
import { Module } from 'vuex';
import { RootState } from '@/store/state';

const defaultState: ExperimentsState = {
  experiments: [],
  dataset: undefined,
  selectedExperimentId: undefined,
  selectedMeta: undefined,
  channels: [],
  selectedAcquisition: undefined,
  selectedMetals: [],
  metalColorMap: {}
};

export const experimentModule: Module<ExperimentsState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
