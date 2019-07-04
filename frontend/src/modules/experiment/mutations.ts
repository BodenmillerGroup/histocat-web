import { IAcquisition, IChannel, IExperiment } from './models';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { RootState } from '@/store/state';

export const mutations = {
  setExperiments(state: ExperimentsState, experiments: IExperiment[]) {
    state.experiments = experiments;
  },
  setTags(state: ExperimentsState, tags: string[]) {
    state.tags = tags;
  },
  setExperiment(state: ExperimentsState, experiment: IExperiment) {
    const items = state.experiments.filter((item: IExperiment) => item.id !== experiment.id);
    items.push(experiment);
    state.experiments = items;
  },
  deleteExperiment(state: ExperimentsState, id: number) {
    state.experiments = state.experiments.filter((item: IExperiment) => item.id !== id);
  },
  setSelectedExperimentId(state: ExperimentsState, id?: number) {
    state.selectedExperimentId = id;
  },
  setChannels(state: ExperimentsState, channels: IChannel[]) {
    state.channels = channels;
  },
  setSelectedAcquisitionId(state: ExperimentsState, id?: number) {
    state.selectedAcquisitionId = id;
  },
  setSelectedMetals(state: ExperimentsState, metals: string[]) {
    state.selectedMetals = metals;
  },
  setMetalColor(state: ExperimentsState, payload: { metal: string, color: string }) {
    state.metalColorMap[payload.metal] = payload.color;
    state.metalColorMap = Object.assign({}, state.metalColorMap);
  },
  setChannelLevels(state: ExperimentsState, payload: { id: number, levels?: { min: number, max: number } }) {
    const channel = state.channels.find((item) => item.id === payload.id);
    if (!channel) {
      return;
    }
    channel.levels = payload.levels;
    state.channels = state.channels.slice();
  },
};

const { commit } = getStoreAccessors<ExperimentsState, RootState>('');

export const commitSetExperiment = commit(mutations.setExperiment);
export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetTags = commit(mutations.setTags);
export const commitDeleteExperiment = commit(mutations.deleteExperiment);
export const commitSetSelectedExperimentId = commit(mutations.setSelectedExperimentId);
export const commitSetChannels = commit(mutations.setChannels);
export const commitSetSelectedAcquisitionId = commit(mutations.setSelectedAcquisitionId);
export const commitSetSelectedMetals = commit(mutations.setSelectedMetals);
export const commitSetMetalColor = commit(mutations.setMetalColor);
export const commitSetChannelLevels = commit(mutations.setChannelLevels);
