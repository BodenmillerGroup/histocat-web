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
import { groupModule } from "@/modules/group";
import { IDatasetUpdate } from "@/modules/datasets/models";
import { uiModule } from "@/modules/ui";

export class DatasetsActions extends Actions<DatasetsState, DatasetsGetters, DatasetsMutations, DatasetsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  ui?: Context<typeof uiModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  projects?: Context<typeof projectsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.ui = uiModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.projects = projectsModule.context(store);
  }

  async getProjectDatasets(projectId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getProjectDatasets(groupId, projectId);
      this.mutations.setEntities(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getDataset(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getDataset(groupId, id);
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
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteDataset(groupId, id);
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

  async uploadDataset(payload: { projectId: number; formData: FormData }) {
    if (!payload.projectId) {
      return;
    }
    try {
      const groupId = this.group?.getters.activeGroupId!;
      await api.uploadDataset(
        this.main!.getters.token,
        groupId,
        payload.projectId,
        payload.formData,
        () => {
          console.log("Upload has started.");
          this.ui!.mutations.setProcessing(true);
        },
        () => {
          console.log("Upload completed successfully.");
          this.ui!.mutations.setProcessing(false);
          this.ui!.mutations.setProcessingProgress(0);
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          this.ui!.mutations.setProcessingProgress(percent);
        },
        () => {}
      );
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
