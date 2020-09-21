import { IDataset } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getProjectDatasets(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/datasets`).json<IDataset[]>();
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
