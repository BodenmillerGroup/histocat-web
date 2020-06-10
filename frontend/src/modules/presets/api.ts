import { IPreset, IPresetCreate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async createPreset(data: IPresetCreate) {
    return ApiManager.api
      .post(`presets`, {
        json: data,
      })
      .json<IPreset>();
  },
  async getExperimentPresets(experimentId: number) {
    return ApiManager.api.get(`experiments/${experimentId}/presets`).json<IPreset[]>();
  },
  async getPreset(id: number) {
    return ApiManager.api.get(`presets/${id}`).json<IPreset>();
  },
  async deletePreset(id: number) {
    return ApiManager.api.delete(`presets/${id}`).json<number>();
  },
};
