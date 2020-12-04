import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { DatasetsState } from ".";
import { api } from "./api";
import { DatasetsGetters } from "./getters";
import { DatasetsMutations } from "./mutations";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_DATASETS, SET_ACTIVE_DATASET_ID } from "./events";
import { groupModule } from "@/modules/group";
import { IDatasetUpdate } from "@/modules/datasets/models";

export class DatasetsActions extends Actions<DatasetsState, DatasetsGetters, DatasetsMutations, DatasetsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  projects?: Context<typeof projectsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.projects = projectsModule.context(store);
  }

  setActiveDatasetId(id: number | null) {
    BroadcastManager.publish(SET_ACTIVE_DATASET_ID, id);
  }

  async getProjectDatasets(projectId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getProjectDatasets(groupId, projectId);
      BroadcastManager.publish(SET_DATASETS, data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getDataset(id: number) {
    try {
      const data = await api.getDataset(id);
      this.mutations.setEntity(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateDataset(payload: { datasetId: number; data: IDatasetUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateDataset(groupId, payload.datasetId, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Dataset successfully updated", color: "success" });
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

  async uploadDataset(payload: { id: number; data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const groupId = this.group?.getters.activeGroupId!;
      await api.uploadDataset(
        this.main!.getters.token,
        groupId,
        payload.id,
        payload.data,
        () => {
          console.log("Upload has started.");
          this.main!.mutations.setProcessing(true);
        },
        () => {
          console.log("Upload completed successfully.");
          this.main!.mutations.setProcessing(false);
          this.main!.mutations.setProcessingProgress(0);
          this.main!.mutations.addNotification({ content: "File successfully uploaded", color: "success" });
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          this.main!.mutations.setProcessingProgress(percent);
        },
        () => {}
      );
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
