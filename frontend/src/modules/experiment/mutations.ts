import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { ExperimentsState } from '.';
import { IExperiment } from './models';

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
  setSelectedAcquisitionId(state: ExperimentsState, id?: number) {
    state.selectedAcquisitionId = id;
  },
  setSelectedMetals(state: ExperimentsState, metals: string[]) {
    state.selectedMetals = metals;
  },
};

const { commit } = getStoreAccessors<ExperimentsState, RootState>('');

export const commitSetExperiment = commit(mutations.setExperiment);
export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetTags = commit(mutations.setTags);
export const commitDeleteExperiment = commit(mutations.deleteExperiment);
export const commitSetSelectedExperimentId = commit(mutations.setSelectedExperimentId);
export const commitSetSelectedAcquisitionId = commit(mutations.setSelectedAcquisitionId);
export const commitSetSelectedMetals = commit(mutations.setSelectedMetals);
