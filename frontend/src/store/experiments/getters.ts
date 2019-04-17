import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { State } from '../state';

export const getters = {
    adminExperiments: (state: ExperimentsState) => state.experiments,
    adminOneExperiment: (state: ExperimentsState) => (experimentId: number) => {
        const filteredExperiments = state.experiments.filter((experiment) => experiment.id === experimentId);
        if (filteredExperiments.length > 0) {
            return { ...filteredExperiments[0] };
        }
    },
};

const { read } = getStoreAccessors<ExperimentsState, State>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readAdminExperiments = read(getters.adminExperiments);
