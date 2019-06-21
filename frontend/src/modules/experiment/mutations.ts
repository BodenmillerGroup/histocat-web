import { IAcquisition, IChannel, IExperiment, IExperimentDataset } from './models';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { RootState } from '@/store/state';

export const mutations = {
  setExperiments(state: ExperimentsState, payload: IExperiment[]) {
    state.experiments = payload;
  },
  setTags(state: ExperimentsState, payload: string[]) {
    state.tags = payload;
  },
  setExperiment(state: ExperimentsState, payload: IExperiment) {
    const items = state.experiments.filter((item: IExperiment) => item.id !== payload.id);
    items.push(payload);
    state.experiments = items;
  },
  deleteExperiment(state: ExperimentsState, id: number) {
    const items = state.experiments.filter((item: IExperiment) => item.id !== id);
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
export const commitSetExperimentDataset = commit(mutations.setExperimentDataset);
export const commitSetSelectedMeta = commit(mutations.setSelectedMeta);
export const commitSetChannels = commit(mutations.setChannels);
export const commitSetSelectedAcquisition = commit(mutations.setSelectedAcquisition);
export const commitSetSelectedMetals = commit(mutations.setSelectedMetals);
export const commitSetMetalColor = commit(mutations.setMetalColor);
export const commitSetChannelLevels = commit(mutations.setChannelLevels);
