import create from "zustand";
import { schema, normalize } from "normalizr";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "../../utils/toaster";
import { IDataset, IDatasetUpdate } from "./models";
import { useGroupsStore } from "../groups";
import { useAuthStore } from "../auth";
import { useMainStore } from "../main";

export const datasetSchema = new schema.Entity("datasets");
export const datasetListSchema = [datasetSchema];

type DatasetsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IDataset };
  activeDatasetId: number | null;

  setActiveDatasetId(id: number | null): void;
  getActiveDataset(): IDataset | null;
  getProjectDatasets(projectId: number): Promise<void>;
  getDataset(id: number): Promise<void>;
  updateDataset(datasetId: number, params: IDatasetUpdate): Promise<void>;
  deleteDataset(id: number): Promise<void>;
  downloadDataset(datasetId: number, filename: string): Promise<void>;
  uploadDataset(id: number, params: any): Promise<void>;
};

export const useDatasetsStore = create<DatasetsState>((set, get) => ({
  ids: [],
  entities: {},
  activeDatasetId: null,

  setActiveDatasetId(id: number | null) {
    set({ activeDatasetId: id });
  },

  getActiveDataset() {
    const activeDatasetId = get().activeDatasetId;
    return activeDatasetId ? get().entities[activeDatasetId] : null;
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

  async updateDataset(datasetId: number, params: IDatasetUpdate) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.updateDataset(groupId, datasetId, params);
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

  async downloadDataset(datasetId: number, filename: string) {
    try {
      const response = await api.downloadDataset(datasetId);
      const blob = await response.blob();
      saveAs(blob, filename);
    } catch (error) {
      displayApiError(error);
    }
  },

  async uploadDataset(id: number, params: any) {
    if (!id) {
      return;
    }
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const token = useAuthStore.getState().token!;
      await api.uploadDataset(
        token,
        groupId,
        id,
        params,
        () => {
          console.log("Upload has started.");
          useMainStore.getState().setProcessing(true);
        },
        () => {
          console.log("Upload completed successfully.");
          useMainStore.getState().setProcessing(false);
          useMainStore.getState().setProcessingProgress(0);
          AppToaster.show({ message: "File successfully uploaded", intent: "success" });
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          useMainStore.getState().setProcessingProgress(percent);
        },
        () => {}
      );
    } catch (error) {
      displayApiError(error);
    }
  },
}));
