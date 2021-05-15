import { IClassifyCellsData, IClassifyCellsSubmission, IRegionChannelData, IRegionStatsSubmission } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async calculateRegionStats(groupId: number, params: IRegionStatsSubmission) {
    return ApiManager.api
      .post(`groups/${groupId}/analysis/region`, {
        json: params,
      })
      .json<IRegionChannelData[]>();
  },
  async classifyCells(groupId: number, params: IClassifyCellsSubmission) {
    return ApiManager.api
      .post(`groups/${groupId}/analysis/classify`, {
        json: params,
        timeout: false,
      })
      .json<IClassifyCellsData>();
  },
};
