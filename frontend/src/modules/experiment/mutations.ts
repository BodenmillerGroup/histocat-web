import { IAcquisition, IChannel, IExperiment, IExperimentDataset } from './models';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '@/store/state';

export const mutations = {
  setExperiments(state: ExperimentsState, payload: IExperiment[]) {
    state.experiments = payload;
  },
  setExperiment(state: ExperimentsState, payload: IExperiment) {
    const items = state.experiments.filter((item: IExperiment) => item.id !== payload.id);
    items.push(payload);
    state.experiments = items;
  },
  setSelectedExperimentId(state: ExperimentsState, payload: { id: number }) {
    state.selectedExperimentId = payload.id;
  },
  setExperimentDataset(state: ExperimentsState, payload: IExperimentDataset) {
    state.dataset = payload;
  },
  setSelectedMeta(state: ExperimentsState, payload: object) {
    state.selectedMeta = payload;
  },
  setChannels(state: ExperimentsState, payload: IChannel[]) {
    state.channels = payload;
  },
  setSelectedAcquisition(state: ExperimentsState, payload: IAcquisition) {
    state.selectedAcquisition = payload;
  },
  setSelectedMetals(state: ExperimentsState, payload: string[]) {
    state.selectedMetals = payload;
  },
};

const { commit } = getStoreAccessors<ExperimentsState, State>('');

export const commitSetExperiment = commit(mutations.setExperiment);
export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetSelectedExperimentId = commit(mutations.setSelectedExperimentId);
export const commitSetExperimentDataset = commit(mutations.setExperimentDataset);
export const commitSetSelectedMeta = commit(mutations.setSelectedMeta);
export const commitSetChannels = commit(mutations.setChannels);
export const commitSetSelectedAcquisition = commit(mutations.setSelectedAcquisition);
export const commitSetSelectedMetals = commit(mutations.setSelectedMetals);
