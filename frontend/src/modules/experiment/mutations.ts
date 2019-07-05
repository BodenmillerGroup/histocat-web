import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { ExperimentsState } from '.';
import { IDataset, IExperiment } from './models';

export const mutations = {
  setExperiments(state: ExperimentsState, experiments: IExperiment[]) {
    state.experiments = experiments;
  },
  setDatasets(state: ExperimentsState, datasets: IDataset[]) {
    state.datasets = datasets;
  },
  setTags(state: ExperimentsState, tags: string[]) {
    state.tags = tags;
  },
  setExperiment(state: ExperimentsState, experiment: IExperiment) {
    const items = state.experiments.filter((item) => item.id !== experiment.id);
    items.push(experiment);
    state.experiments = items;
  },
  setDataset(state: ExperimentsState, dataset: IDataset) {
    const items = state.datasets.filter((item) => item.id !== dataset.id);
    items.push(dataset);
    state.datasets = items;
  },
  deleteExperiment(state: ExperimentsState, id: number) {
    state.experiments = state.experiments.filter((item) => item.id !== id);
  },
  deleteDataset(state: ExperimentsState, id: number) {
    state.datasets = state.datasets.filter((item) => item.id !== id);
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

export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetExperiment = commit(mutations.setExperiment);
export const commitDeleteExperiment = commit(mutations.deleteExperiment);

export const commitSetDatasets = commit(mutations.setDatasets);
export const commitSetDataset = commit(mutations.setDataset);
export const commitDeleteDataset = commit(mutations.deleteDataset);

export const commitSetTags = commit(mutations.setTags);
export const commitSetSelectedExperimentId = commit(mutations.setSelectedExperimentId);
export const commitSetSelectedAcquisitionId = commit(mutations.setSelectedAcquisitionId);
export const commitSetSelectedMetals = commit(mutations.setSelectedMetals);
