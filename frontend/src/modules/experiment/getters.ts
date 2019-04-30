import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '@/store/state';

export const getters = {
  adminExperiments: (state: ExperimentsState) => state.experiments,
  adminOneExperiment: (state: ExperimentsState) => (id: number) => {
    return state.experiments.find((item) => item.id === id);
  },
  activeExperiment: (state: ExperimentsState) => {
    return state.experiments.find((item) => item.id === state.activeExperimentId);
  },
  dataset: (state: ExperimentsState) => {
    return state.dataset;
  },
  activeMeta: (state: ExperimentsState) => {
    return state.activeMeta;
  },
};

const {read} = getStoreAccessors<ExperimentsState, State>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readAdminExperiments = read(getters.adminExperiments);
export const readActiveExperiment = read(getters.activeExperiment);
export const readExperimentDataset = read(getters.dataset);
export const readActiveMeta = read(getters.activeMeta);
