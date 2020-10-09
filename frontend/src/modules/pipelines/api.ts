import { IPipeline, IPipelineCreate, IPipelineUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createPipeline(groupId: number, data: IPipelineCreate) {
    return ApiManager.api
      .post(`groups/${groupId}/pipelines`, {
        json: data,
      })
      .json<IPipeline>();
  },
  async updatePipeline(groupId: number, pipelineId: number, data: IPipelineUpdate) {
    return ApiManager.api
      .put(`groups/${groupId}/pipelines/${pipelineId}`, {
        json: data,
      })
      .json<IPipeline>();
  },
  async getProjectPipelines(projectId: number) {
    return ApiManager.api.get(`projects/${projectId}/pipelines`).json<IPipeline[]>();
  },
  async getPipeline(id: number) {
    return ApiManager.api.get(`pipelines/${id}`).json<IPipeline>();
  },
  async deletePipeline(id: number) {
    return ApiManager.api.delete(`pipelines/${id}`).json<number>();
  },
};
