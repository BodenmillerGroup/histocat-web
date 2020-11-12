import { IRawResultData, IResult, IResultUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getDatasetResults(groupId: number, datasetId: number) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${datasetId}/results`).json<IResult[]>();
  },
  async getResult(groupId: number, resultId: number) {
    return ApiManager.api.get(`groups/${groupId}/results/${resultId}`).json<IResult>();
  },
  async getResultData(groupId: number, resultId: number, colorsType?: string, colorsName?: string) {
    let url = `groups/${groupId}/results/${resultId}/data`;
    if (colorsType && colorsName) {
      url += `?colors_type=${colorsType}&colors_name=${colorsName}`;
    }
    return ApiManager.api.get(url).json<IRawResultData>();
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
