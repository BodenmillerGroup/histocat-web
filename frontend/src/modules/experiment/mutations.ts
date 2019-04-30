import { IExperiment, IExperimentDataset } from './models';
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
  setActiveExperimentId(state: ExperimentsState, payload: { id: number }) {
    state.activeExperimentId = payload.id;
  },
  setExperimentDataset(state: ExperimentsState, payload: IExperimentDataset) {
    state.dataset = payload;
  },
  setActiveMeta(state: ExperimentsState, payload: object) {
    state.activeMeta = payload;
  },
};

const {commit} = getStoreAccessors<ExperimentsState, State>('');

export const commitSetExperiment = commit(mutations.setExperiment);
export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetActiveExperimentId = commit(mutations.setActiveExperimentId);
export const commitSetExperimentDataset = commit(mutations.setExperimentDataset);
export const commitSetActiveMeta = commit(mutations.setActiveMeta);
