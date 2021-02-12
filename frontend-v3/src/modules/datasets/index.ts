import create from "zustand";
import { schema, normalize } from "normalizr";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { isEqual } from "lodash-es";
import { AppToaster } from "../../utils/toaster";
import { IDataset, IDatasetUpdate } from "./models";
import { useGroupsStore } from "../groups";

export const datasetSchema = new schema.Entity("datasets");
export const datasetListSchema = [datasetSchema];

type DatasetsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IDataset };
  activeDatasetId: number | null;

  setActiveDatasetId(id: number | null): void;
};

export const useDatasetsStore = create<DatasetsState>((set, get) => ({
  ids: [],
  entities: {},
  activeDatasetId: null,

  setActiveDatasetId(id: number | null) {
    set({ activeDatasetId: id });
  },

  async getProjectDatasets(projectId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getProjectDatasets(groupId, projectId);
      if (data) {
        const normalizedData = normalize<IDataset>(data, datasetListSchema);
        set({
          ids: normalizedData.result,
          entities: normalizedData.entities.datasets ? normalizedData.entities.datasets : {},
        });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getDataset(id: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getDataset(groupId, id);
      if (data) {
        const existingId = get().ids.find((id) => id === data.id);
        if (!existingId) {
          set({ ids: get().ids.concat(data.id) });
        }
        set({ entities: { ...get().entities, [data.id]: data } });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateDataset(payload: { datasetId: number; data: IDatasetUpdate }) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.updateDataset(groupId, payload.datasetId, payload.data);
      if (data) {
        set({ entities: { ...get().entities, [data.id]: data } });
        AppToaster.show({ message: "Dataset successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteDataset(id: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.deleteDataset(groupId, id);
      if (data) {
        const entities = { ...get().entities };
        delete entities[id];
        set({ ids: get().ids.filter((item) => item !== id), entities: entities });
        AppToaster.show({ message: "Dataset successfully deleted", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async downloadDataset(payload: { datasetId: number; filename: string }) {
    try {
      const response = await api.downloadDataset(payload.datasetId);
      const blob = await response.blob();
      saveAs(blob, payload.filename);
    } catch (error) {
      displayApiError(error);
    }
  },

  async uploadDataset(payload: { id: number; data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
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
      displayApiError(error);
    }
  },
}));
