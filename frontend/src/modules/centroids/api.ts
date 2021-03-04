import { ICentroidsData, ICentroidsSubmission } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getCentroids(groupId: number, params: ICentroidsSubmission) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${params.datasetId}/centroids`).json<ICentroidsData>();
  },
};
