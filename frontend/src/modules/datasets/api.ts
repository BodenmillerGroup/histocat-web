import { IDataset, IDatasetCreate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createDataset(data: IDatasetCreate) {
    return ApiManager.api
      .post(`datasets`, {
        json: data,
      })
      .json<IDataset>();
  },
  async uploadDataset(experimentId: number, data) {
    return ApiManager.api.post(`datasets/experiment/${experimentId}/upload`, {
      body: data,
      timeout: false,
    });
  },
  async getExperimentDatasets(groupId: number, experimentId: number) {
    return ApiManager.api.get(`groups/${groupId}/experiments/${experimentId}/datasets`).json<IDataset[]>();
  },
  async deleteDataset(id: number) {
    return ApiManager.api.delete(`datasets/${id}`).json();
  },
  async downloadDataset(id: number) {
    return ApiManager.api.get(`datasets/${id}/download`, {
      timeout: false,
    });
  },
};
