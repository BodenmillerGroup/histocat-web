import { IRegionChannelData, IRegionStatsSubmission } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async calculateRegionStats(groupId: number, params: IRegionStatsSubmission) {
    return ApiManager.api
      .post(`groups/${groupId}/analysis/region/stats`, {
        json: params,
      })
      .json<IRegionChannelData[]>();
  },
};
