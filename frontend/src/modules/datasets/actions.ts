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
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_DATASETS } from "@/modules/datasets/events";

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
      const data = await api.getExperimentDatasets(experimentId);
      BroadcastManager.publish(SET_DATASETS, data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteDataset(id: number) {
    try {
      const data = await api.deleteDataset(id);
      this.mutations.deleteEntity(id);
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
        for (const acquisition of slide.acquisitions) {
          if (acquisitionIds.includes(acquisition.id)) {
            for (const channel of Object.values(acquisition.channels)) {
              if (metals.includes(channel.name)) {
                const settings = this.settings!.getters.getChannelSettings(acquisition.id, channel.name);
                if (settings) {
                  channelsSettings.push(settings);
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
      const data = await api.createDataset(params);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Dataset successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async downloadDataset(payload: { datasetId: number; filename: string }) {
    try {
      const response = await api.downloadDataset(payload.datasetId);
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
      await api.uploadDataset(payload.experimentId, payload.data);
      this.main!.mutations.addNotification({ content: "Dataset successfully uploaded", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
