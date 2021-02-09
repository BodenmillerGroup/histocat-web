import { IPreset, IPresetCreate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createPreset(groupId: number, data: IPresetCreate) {
    return ApiManager.api
      .post(`groups/${groupId}/presets`, {
        json: data,
      })
      .json<IPreset>();
  },
  async getProjectPresets(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/presets`).json<IPreset[]>();
  },
  async getPreset(groupId: number, id: number) {
    return ApiManager.api.get(`groups/${groupId}/presets/${id}`).json<IPreset>();
  },
  async deletePreset(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/presets/${id}`).json<number>();
  },
};
