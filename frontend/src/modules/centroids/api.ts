import {
  ICentroidsData,
  ICentroidsSubmission,
} from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getCentroids(params: ICentroidsSubmission) {
    return ApiManager.api
      .get(`datasets/${params.datasetId}/centroids`)
      .json<ICentroidsData>();
  },
};
