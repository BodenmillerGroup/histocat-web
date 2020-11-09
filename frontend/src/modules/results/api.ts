import { IResult, IResultUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getDatasetResults(groupId: number, datasetId: number) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${datasetId}/results`).json<IResult[]>();
  },
  async getResult(groupId: number, resultId: number) {
    return ApiManager.api.get(`groups/${groupId}/results/${resultId}`).json<IResult>();
  },
  async updateResult(groupId: number, resultId: number, data: IResultUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/results/${resultId}`, {
        json: data,
      })
      .json<IResult>();
  },
  async deleteResult(groupId: number, resultId: number) {
    return ApiManager.api.delete(`groups/${groupId}/results/${resultId}`).json();
  },
};
