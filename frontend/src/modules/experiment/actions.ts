import { mainModule } from '@/modules/main';
import { settingsModule } from '@/modules/settings';
import { IChannelSettings } from '@/modules/settings/models';
import { Store } from 'vuex';
import { Actions, Context } from 'vuex-smart-module';
import { ExperimentState } from '.';
import { api } from './api';
import { ExperimentGetters } from './getters';
import { IDatasetCreate, IExperimentCreate, IExperimentUpdate } from './models';
import { ExperimentMutations } from './mutations';

export class ExperimentActions extends Actions<ExperimentState, ExperimentGetters, ExperimentMutations, ExperimentActions> {

  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
  }

  async getExperiments() {
    try {
      const data = await api.getExperiments(this.main!.getters.token);
      if (data) {
        this.mutations.setExperiments(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getTags() {
    try {
      const data = await api.getTags(this.main!.getters.token);
      if (data) {
        this.mutations.setTags(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateExperiment(payload: { id: number, data: IExperimentUpdate }) {
    try {
      const notification = { content: 'saving', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.updateExperiment(this.main!.getters.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.setExperiment(data as any);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Experiment successfully updated', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteExperiment(id: number) {
    try {
      const notification = { content: 'deleting', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteExperiment(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.deleteExperiment(id);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Experiment successfully deleted', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createExperiment(payload: IExperimentCreate) {
    try {
      const notification = { content: 'saving', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createExperiment(this.main!.getters.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.setExperiment(data as any);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Experiment successfully created', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async uploadSlide(payload: { id: number, data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const notification = { content: 'saving', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const response = (await Promise.all([
        api.uploadSlide(this.main!.getters.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'File successfully uploaded', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getExperimentData(id: number) {
    try {
      const data = await api.getExperimentData(this.main!.getters.token, id);
      if (data) {
        this.mutations.setExperiment(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStats(id: number) {
    try {
      return await api.getChannelStats(this.main!.getters.token, id);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getOwnDatasets(experimentId: number) {
    try {
      const data = await api.getOwnDatasets(this.main!.getters.token, experimentId);
      this.mutations.setDatasets(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteDataset(id: number) {
    try {
      const notification = { content: 'deleting', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteDataset(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.deleteDataset(id);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Dataset successfully deleted', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createDataset(payload: { name: string, description: string }) {
    const experimentId = this.getters.activeExperimentId;
    const acquisitionIds = this.getters.selectedAcquisitionIds;
    if (acquisitionIds.length === 0) {
      this.main!.mutations.addNotification({ content: 'Please select at least one acquisition', color: 'warning' });
      return;
    }
    const metals = this.getters.selectedMetals;
    if (metals.length === 0) {
      this.main!.mutations.addNotification({ content: 'Please select at least one channel', color: 'warning' });
      return;
    }
    const channelsSettings: IChannelSettings[] = [];

    const experiment = this.getters.activeExperiment;
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            for (const acquisition of roi.acquisitions) {
              if (acquisitionIds.includes(acquisition.id)) {
                for (const channel of acquisition.channels) {
                  if (metals.includes(channel.metal)) {
                    const settings = this.settings!.getters.channelSettings(channel.id);
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
      input: {
        acquisition_ids: acquisitionIds,
        metals: metals,
        channel_settings: channelsSettings,
      },
    };

    try {
      const notification = { content: 'saving', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createDataset(this.main!.getters.token, params),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.setDataset(data as any);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Dataset successfully created', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
