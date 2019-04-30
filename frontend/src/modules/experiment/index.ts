import { mutations } from './mutations';
import { getters } from './getters';
import { actions } from './actions';
import { ExperimentsState } from './state';

const defaultState: ExperimentsState = {
  experiments: [],
  dataset: undefined,
  activeExperimentId: undefined,
  activeMeta: undefined,
};

export const experimentModule = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
