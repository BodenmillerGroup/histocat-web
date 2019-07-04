import { api } from './api';
import { ActionContext } from 'vuex';
import { IExperimentCreate, IExperimentUpdate } from './models';
import { RootState } from '@/store/state';
import { ExperimentsState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import {
  commitDeleteExperiment,
  commitSetExperiment,
  commitSetExperiments,
  commitSetTags,
} from './mutations';
import { dispatchCheckApiError } from '@/modules/main/actions';
import { commitAddNotification, commitRemoveNotification } from '@/modules/main/mutations';

type ExperimentContext = ActionContext<ExperimentsState, RootState>;

export const actions = {
  async actionGetExperiments(context: ExperimentContext) {
    try {
      const data = await api.getExperiments(context.rootState.main.token);
      if (data) {
        commitSetExperiments(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetTags(context: ExperimentContext) {
    try {
      const data = await api.getTags(context.rootState.main.token);
      if (data) {
        commitSetTags(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionUpdateExperiment(context: ExperimentContext, payload: { id: number, data: IExperimentUpdate }) {
    try {
      const notification = { content: 'saving', showProgress: true };
      commitAddNotification(context, notification);
      const data = (await Promise.all([
        api.updateExperiment(context.rootState.main.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetExperiment(context, data as any);
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'Experiment successfully updated', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionDeleteExperiment(context: ExperimentContext, id: number) {
    try {
      const notification = { content: 'deleting', showProgress: true };
      commitAddNotification(context, notification);
      const data = (await Promise.all([
        api.deleteExperiment(context.rootState.main.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitDeleteExperiment(context, id);
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'Experiment successfully deleted', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionCreateExperiment(context: ExperimentContext, payload: IExperimentCreate) {
    try {
      const notification = { content: 'saving', showProgress: true };
      commitAddNotification(context, notification);
      const data = (await Promise.all([
        api.createExperiment(context.rootState.main.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetExperiment(context, data as any);
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'Experiment successfully created', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionUploadSlide(context: ExperimentContext, payload: { id: number, data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const notification = { content: 'saving', showProgress: true };
      commitAddNotification(context, notification);
      const response = (await Promise.all([
        api.uploadSlide(context.rootState.main.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'File successfully uploaded', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetExperimentDataset(context: ExperimentContext, id: number) {
    try {
      const data = await api.getExperimentDataset(context.rootState.main.token, id);
      if (data) {
        commitSetExperiment(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetChannelImage(context: ExperimentContext, id: number) {
    try {
      const data = await api.getChannelImage(context.rootState.main.token, id);
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetChannelStats(context: ExperimentContext, id: number) {
    try {
      const data = await api.getChannelStats(context.rootState.main.token, id);
      return data;
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
};

const { dispatch } = getStoreAccessors<ExperimentsState, RootState>('');

export const dispatchCreateExperiment = dispatch(actions.actionCreateExperiment);
export const dispatchGetExperiments = dispatch(actions.actionGetExperiments);
export const dispatchGetTags = dispatch(actions.actionGetTags);
export const dispatchUpdateExperiment = dispatch(actions.actionUpdateExperiment);
export const dispatchDeleteExperiment = dispatch(actions.actionDeleteExperiment);
export const dispatchUploadSlide = dispatch(actions.actionUploadSlide);
export const dispatchGetExperimentDataset = dispatch(actions.actionGetExperimentDataset);
export const dispatchGetChannelImage = dispatch(actions.actionGetChannelImage);
export const dispatchGetChannelStats = dispatch(actions.actionGetChannelStats);
