import { IDatasetCreate } from "@/modules/datasets/models";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { IChannelSettings } from "@/modules/settings/models";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { DatasetState } from ".";
import { api } from "./api";
import { DatasetGetters } from "./getters";
import { DatasetMutations } from "./mutations";

export class DatasetActions extends Actions<DatasetState, DatasetGetters, DatasetMutations, DatasetActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;
  experiment?: Context<typeof experimentModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
    this.experiment = experimentModule.context(store);
  }

  async getExperimentDatasets(experimentId: number) {
    try {
      const data = await api.getExperimentDatasets(this.main!.getters.token, experimentId);
      this.mutations.setDatasets(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteDataset(id: number) {
    try {
      const notification = { content: "deleting", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteDataset(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.deleteDataset(id);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Dataset successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createDataset(payload: { name: string; description: string }) {
    const experimentId = this.experiment!.getters.activeExperimentId;
    const acquisitionIds = this.experiment!.getters.selectedAcquisitionIds;
    if (acquisitionIds.length === 0) {
      this.main!.mutations.addNotification({ content: "Please select at least one acquisition", color: "warning" });
      return;
    }
    const metals = this.experiment!.getters.selectedMetals;
    if (metals.length === 0) {
      this.main!.mutations.addNotification({ content: "Please select at least one channel", color: "warning" });
      return;
    }
    const channelsSettings: IChannelSettings[] = [];

    const experiment = this.experiment!.getters.activeExperiment;
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            for (const acquisition of roi.acquisitions) {
              if (acquisitionIds.includes(acquisition.id)) {
                for (const channel of acquisition.channels) {
                  if (metals.includes(channel.metal)) {
                    const settings = this.settings!.getters.getChannelSettings(channel.id);
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
        channel_settings: channelsSettings
      }
    };

    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createDataset(this.main!.getters.token, params),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setDataset(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Dataset successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async downloadDataset(payload: { datasetId: number; filename: string }) {
    try {
      const response = await api.downloadDataset(this.main!.getters.token, payload.datasetId);
      const blob = await response.blob();
      saveAs(blob, payload.filename);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async uploadDataset(payload: { experimentId?: number; data: any }) {
    if (!payload.experimentId) {
      return;
    }
    try {
      const notification = { content: "uploading", showProgress: true };
      this.main!.mutations.addNotification(notification);
      await api.uploadDataset(this.main!.getters.token, payload.experimentId, payload.data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Dataset successfully uploaded", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
