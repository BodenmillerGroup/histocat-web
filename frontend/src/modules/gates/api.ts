import { IGate, IGateCreate, IGateUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createGate(data: IGateCreate) {
    return ApiManager.api
      .post(`gates`, {
        json: data,
      })
      .json<IGate>();
  },
  async updateGate(groupId: number, gateId: number, data: IGateUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/gates/${gateId}`, {
        json: data,
      })
      .json<IGate>();
  },
  async getDatasetGates(datasetId: number) {
    return ApiManager.api.get(`datasets/${datasetId}/gates`).json<IGate[]>();
  },
  async getGate(id: number) {
    return ApiManager.api.get(`gates/${id}`).json<IGate>();
  },
  async deleteGate(id: number) {
    return ApiManager.api.delete(`gates/${id}`).json<number>();
  },
};
