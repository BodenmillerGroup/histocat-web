import { IExperiment } from './models';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '@/store/state';

export const mutations = {
  setExperiments(state: ExperimentsState, payload: IExperiment[]) {
    state.experiments = payload;
  },
  setExperiment(state: ExperimentsState, payload: IExperiment) {
    const experiments = state.experiments.filter((experiment: IExperiment) => experiment.id !== payload.id);
    experiments.push(payload);
    state.experiments = experiments;
  },
  setActiveExperimentId(state: ExperimentsState, payload: { id: number }) {
    state.activeExperimentId = payload.id;
  },
};

const {commit} = getStoreAccessors<ExperimentsState, State>('');

export const commitSetExperiment = commit(mutations.setExperiment);
export const commitSetExperiments = commit(mutations.setExperiments);
export const commitSetActiveExperimentId = commit(mutations.setActiveExperimentId);
