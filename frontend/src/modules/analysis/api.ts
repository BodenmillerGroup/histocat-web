import { IRegionChannelData, IRegionStatsSubmission } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async calculateRegionStats(params: IRegionStatsSubmission) {
    return ApiManager.api
      .post(`analysis/region/stats`, {
        json: params,
      })
      .json<IRegionChannelData[]>();
  },
};
