import { IDataset, IDatasetUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getProjectDatasets(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/datasets`).json<IDataset[]>();
  },
  async getDataset(datasetId: number) {
    return ApiManager.api.get(`datasets/${datasetId}`).json<IDataset>();
  },
  async updateDataset(groupId: number, datasetId: number, data: IDatasetUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/datasets/${datasetId}`, {
        json: data,
      })
      .json<IDataset>();
  },
  async deleteDataset(datasetId: number) {
    return ApiManager.api.delete(`datasets/${datasetId}`).json();
  },
  async downloadDataset(datasetId: number) {
    return ApiManager.api.get(`datasets/${datasetId}/download`, {
      timeout: false,
    });
  },
};
