import { apiUrl } from "@/env";
import ky from "ky";
import { IDataset, IDatasetCreate } from "./models";

export const api = {
  async createDataset(token: string, data: IDatasetCreate) {
    return ky
      .post(`${apiUrl}/api/v1/datasets/`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IDataset>();
  },
  async uploadDataset(token: string, experimentId: number, data) {
    return ky.post(`${apiUrl}/api/v1/datasets/experiment/${experimentId}/upload`, {
      body: data,
      headers: {
        Authorization: `Bearer ${token}`
      },
      timeout: false
    });
  },
  async getExperimentDatasets(token: string, experimentId: number) {
    return ky
      .get(`${apiUrl}/api/v1/datasets/experiment/${experimentId}`, {
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IDataset[]>();
  },
  async deleteDataset(token: string, id: number) {
    return ky
      .delete(`${apiUrl}/api/v1/datasets/${id}`, {
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json();
  },
  async downloadDataset(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/datasets/${id}/download`, {
      headers: {
        Authorization: `Bearer ${token}`
      },
      timeout: false
    });
  }
};
