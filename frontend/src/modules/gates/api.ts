import { IGate, IGateCreate, IGateUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createGate(groupId: number, data: IGateCreate) {
    return ApiManager.api
      .post(`groups/${groupId}/gates`, {
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
  async getDatasetGates(groupId: number, datasetId: number) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${datasetId}/gates`).json<IGate[]>();
  },
  async getGate(groupId: number, id: number) {
    return ApiManager.api.get(`groups/${groupId}/gates/${id}`).json<IGate>();
  },
  async deleteGate(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/gates/${id}`).json<number>();
  },
};
