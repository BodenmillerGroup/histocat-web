import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { DatasetState } from ".";
import { api } from "./api";
import { DatasetGetters } from "./getters";
import { DatasetMutations } from "./mutations";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_DATASETS, SET_ACTIVE_DATASET_ID } from "./events";
import { groupModule } from "@/modules/group";

export class DatasetActions extends Actions<DatasetState, DatasetGetters, DatasetMutations, DatasetActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  experiment?: Context<typeof experimentModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.experiment = experimentModule.context(store);
  }

  setActiveDatasetId(id: number | null) {
    BroadcastManager.publish(SET_ACTIVE_DATASET_ID, id);
  }

  async getExperimentDatasets(experimentId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getExperimentDatasets(groupId, experimentId);
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
