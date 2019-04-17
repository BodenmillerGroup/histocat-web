import { api } from './api';
import { ActionContext } from 'vuex';
import { IExperimentCreate, IExperimentUpdate } from './models';
import { State } from '@/store/state';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { commitSetExperiment, commitSetExperiments } from './mutations';
import { dispatchCheckApiError } from '@/modules/main/actions';
import { commitAddNotification, commitRemoveNotification } from '@/modules/main/mutations';

type MainContext = ActionContext<ExperimentsState, State>;

export const actions = {
  async actionGetExperiments(context: MainContext) {
    try {
      const response = await api.getExperiments(context.rootState.main.token);
      if (response) {
        commitSetExperiments(context, response.data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionUpdateExperiment(context: MainContext, payload: { id: number, data: IExperimentUpdate }) {
    try {
      const loadingNotification = {content: 'saving', showProgress: true};
      commitAddNotification(context, loadingNotification);
      const response = (await Promise.all([
        api.updateExperiment(context.rootState.main.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetExperiment(context, response.data);
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, {content: 'Experiment successfully updated', color: 'success'});
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionCreateExperiment(context: MainContext, payload: IExperimentCreate) {
    try {
      const loadingNotification = {content: 'saving', showProgress: true};
      commitAddNotification(context, loadingNotification);
      const response = (await Promise.all([
        api.createExperiment(context.rootState.main.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetExperiment(context, response.data);
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, {content: 'Experiment successfully created', color: 'success'});
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
};

const {dispatch} = getStoreAccessors<ExperimentsState, State>('');

export const dispatchCreateExperiment = dispatch(actions.actionCreateExperiment);
export const dispatchGetExperiments = dispatch(actions.actionGetExperiments);
export const dispatchUpdateExperiment = dispatch(actions.actionUpdateExperiment);
