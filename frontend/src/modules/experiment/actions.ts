import { api } from './api';
import { ActionContext } from 'vuex';
import { IExperimentCreate, IExperimentUpdate } from './models';
import { State } from '@/store/state';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { commitSetActiveExperimentId, commitSetExperiment, commitSetExperiments } from './mutations';
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
  async actionUploadSlide(context: MainContext, payload: { id: number | undefined, data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const loadingNotification = {content: 'saving', showProgress: true};
      commitAddNotification(context, loadingNotification);
      const response = (await Promise.all([
        api.uploadSlide(context.rootState.main.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, {content: 'File successfully uploaded', color: 'success'});
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetExperiment(context: MainContext, payload: { id: number }) {
    try {
      const response = await api.getExperiment(context.rootState.main.token, payload.id);
      if (response) {
        commitSetExperiment(context, response.data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionSetActiveExperimentId(context: MainContext, payload: { id: number }) {
    commitSetActiveExperimentId(context, payload);
  },
};

const {dispatch} = getStoreAccessors<ExperimentsState, State>('');

export const dispatchCreateExperiment = dispatch(actions.actionCreateExperiment);
export const dispatchGetExperiments = dispatch(actions.actionGetExperiments);
export const dispatchUpdateExperiment = dispatch(actions.actionUpdateExperiment);
export const dispatchUploadSlide = dispatch(actions.actionUploadSlide);
export const dispatchGetExperiment = dispatch(actions.actionGetExperiment);
export const dispatchSetActiveExperimentId = dispatch(actions.actionSetActiveExperimentId);
