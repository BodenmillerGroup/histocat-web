import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '@/store/state';

export const getters = {
  adminExperiments: (state: ExperimentsState) => state.experiments,
  adminOneExperiment: (state: ExperimentsState) => (id: number) => {
    return state.experiments.find((item) => item.id === id);
  },
  selectedExperiment: (state: ExperimentsState) => {
    return state.experiments.find((item) => item.id === state.selectedExperimentId);
  },
  dataset: (state: ExperimentsState) => {
    return state.dataset;
  },
  selectedMeta: (state: ExperimentsState) => {
    return state.selectedMeta;
  },
  channels: (state: ExperimentsState) => {
    return state.channels;
  },
  selectedAcquisition: (state: ExperimentsState) => {
    return state.selectedAcquisition;
  },
  selectedMetals: (state: ExperimentsState) => {
    return state.selectedMetals;
  },
  selectedChannels: (state: ExperimentsState) => {
    const channels = state.channels.filter((channel) => {
      if (state.selectedMetals.includes(channel.metal)) {
        return channel;
      }
    });
    return channels;
  },
};

const { read } = getStoreAccessors<ExperimentsState, State>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readAdminExperiments = read(getters.adminExperiments);
export const readSelectedExperiment = read(getters.selectedExperiment);
export const readExperimentDataset = read(getters.dataset);
export const readSelectedMeta = read(getters.selectedMeta);
export const readChannels = read(getters.channels);
export const readSelectedAcquisition = read(getters.selectedAcquisition);
export const readSelectedMetals = read(getters.selectedMetals);
export const readSelectedChannels = read(getters.selectedChannels);
