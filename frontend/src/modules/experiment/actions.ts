import {
  readActiveExperiment,
  readActiveExperimentId,
  readSelectedAcquisitionIds,
  readSelectedMetals,
} from '@/modules/experiment/getters';
import { dispatchCheckApiError } from '@/modules/main/actions';
import { commitAddNotification, commitRemoveNotification } from '@/modules/main/mutations';
import { IChannelSettings } from '@/modules/settings/models';
import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { ActionContext } from 'vuex';
import { ExperimentsState } from '.';
import { api } from './api';
import { IDatasetCreate, IExperimentCreate, IExperimentUpdate } from './models';
import {
  commitDeleteDataset,
  commitDeleteExperiment,
  commitSetDataset,
  commitSetDatasets,
  commitSetExperiment,
  commitSetExperiments,
  commitSetTags,
} from './mutations';

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
  async actionGetExperimentData(context: ExperimentContext, id: number) {
    try {
      const data = await api.getExperimentData(context.rootState.main.token, id);
      if (data) {
        commitSetExperiment(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetChannelStats(context: ExperimentContext, id: number) {
    try {
      return await api.getChannelStats(context.rootState.main.token, id);
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionGetOwnDatasets(context: ExperimentContext, experimentId: number) {
    try {
      const data = await api.getOwnDatasets(context.rootState.main.token, experimentId);
      commitSetDatasets(context, data);
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionDeleteDataset(context: ExperimentContext, id: number) {
    try {
      const notification = { content: 'deleting', showProgress: true };
      commitAddNotification(context, notification);
      const data = (await Promise.all([
        api.deleteDataset(context.rootState.main.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitDeleteDataset(context, id);
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'Dataset successfully deleted', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionCreateDataset(context: ExperimentContext, payload: { name: string, description: string }) {
    const experimentId = readActiveExperimentId(context);
    const acquisitionIds = readSelectedAcquisitionIds(context);
    const metals = readSelectedMetals(context);
    const channelsSettings: IChannelSettings[] = [];

    const experiment = readActiveExperiment(context);
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            for (const acquisition of roi.acquisitions) {
              if (acquisitionIds.includes(acquisition.id)) {
                for (const channel of acquisition.channels) {
                  if (metals.includes(channel.metal)) {
                    const settings = context.getters.channelSettings(channel.id);
                    if (settings) {
                      channelsSettings.push(settings);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    const params: IDatasetCreate = {
      experiment_id: experimentId!,
      name: payload.name,
      description: payload.description,
      meta: {
        input: {
          acquisition_ids: acquisitionIds,
          metals: metals,
          channel_settings: channelsSettings,
        },
      },
    };

    try {
      const notification = { content: 'saving', showProgress: true };
      commitAddNotification(context, notification);
      const data = (await Promise.all([
        api.createDataset(context.rootState.main.token, params),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetDataset(context, data as any);
      commitRemoveNotification(context, notification);
      commitAddNotification(context, { content: 'Dataset successfully created', color: 'success' });
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
export const dispatchGetExperimentData = dispatch(actions.actionGetExperimentData);
export const dispatchGetChannelStats = dispatch(actions.actionGetChannelStats);
export const dispatchCreateDataset = dispatch(actions.actionCreateDataset);
export const dispatchGetOwnDatasets = dispatch(actions.actionGetOwnDatasets);
export const dispatchDeleteDataset = dispatch(actions.actionDeleteDataset);
