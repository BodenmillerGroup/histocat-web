import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '@/store/state';

export const getters = {
  adminExperiments: (state: ExperimentsState) => state.experiments,
  adminOneExperiment: (state: ExperimentsState) => (experimentId: number) => {
    const filteredExperiments = state.experiments.filter((experiment) => experiment.id === experimentId);
    if (filteredExperiments.length > 0) {
      return {...filteredExperiments[0]};
    }
  },
  activeExperiment: (state: ExperimentsState) => {
    const activeExperiment = state.experiments.find((experiment) => experiment.id === state.activeExperimentId);
    return activeExperiment;
  },
};

const {read} = getStoreAccessors<ExperimentsState, State>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readAdminExperiments = read(getters.adminExperiments);
export const readActiveExperiment = read(getters.activeExperiment);
