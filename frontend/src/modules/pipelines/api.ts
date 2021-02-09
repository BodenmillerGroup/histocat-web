import { IPipeline, IPipelineCreate, IPipelineUpdate, IProcessPipeline } from "./models";
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
      .patch(`groups/${groupId}/pipelines/${pipelineId}`, {
        json: data,
      })
      .json<IPipeline>();
  },
  async getProjectPipelines(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/pipelines`).json<IPipeline[]>();
  },
  async getPipeline(groupId: number, id: number) {
    return ApiManager.api.get(`groups/${groupId}/pipelines/${id}`).json<IPipeline>();
  },
  async deletePipeline(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/pipelines/${id}`).json<number>();
  },
  async processPipeline(groupId: number, data: IProcessPipeline) {
    return ApiManager.api
      .post(`groups/${groupId}/pipelines/process`, {
        json: data,
      })
      .json();
  },
};
