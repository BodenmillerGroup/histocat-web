import { ApiManager } from "@/utils/api";
import { ISegmentationSubmission } from "./models";

export const api = {
  async processSegmentation(groupId: number, projectId: number, params: ISegmentationSubmission) {
    return ApiManager.api
      .post(`groups/${groupId}/projects/${projectId}/segmentation/process`, {
        json: params,
      })
      .json<boolean>();
  },
};
