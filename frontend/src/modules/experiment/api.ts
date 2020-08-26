import { apiUrl } from "@/env";
import {
  IChannelStack,
  IChannelStats,
  IExperiment,
  IExperimentCreate,
  IExperimentData,
  IExperimentUpdate,
  ISlide,
} from "./models";
import { ApiManager } from "@/utils/api";

const cacheAvailable = false; // 'caches' in self;

export const api = {
  async getGroupExperiments(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/experiments`).json<IExperiment[]>();
  },
  async getExperimentTags(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/tags`).json<string[]>();
  },
  async updateExperiment(groupId: number, experimentId: number, data: IExperimentUpdate) {
    return ApiManager.api
      .put(`groups/${groupId}/experiments/${experimentId}`, {
        json: data,
      })
      .json<IExperiment>();
  },
  async createExperiment(data: IExperimentCreate) {
    return ApiManager.api
      .post(`groups/${data.group_id}/experiments`, {
        json: data,
      })
      .json<IExperiment>();
  },
  async upload(
    token: string,
    groupId: number,
    experimentId: number,
    data,
    onLoadStart: () => void,
    onLoad: () => void,
    onProgress: (event: ProgressEvent) => void,
    onError: () => void
  ) {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", `${apiUrl}/groups/${groupId}/experiments/${experimentId}/upload`);
    xhr.setRequestHeader("Authorization", `Bearer ${token}`);

    xhr.upload.onloadstart = onLoadStart;
    xhr.upload.onprogress = onProgress;
    xhr.upload.onload = onLoad;
    xhr.upload.onerror = function () {
      console.log(`Error during file upload: ${xhr.status}.`);
      onError();
    };

    xhr.send(data);

    // return ky.post(`${apiUrl}/api/v1/experiments/${id}/upload`, {
    //   body: data,
    //   headers: {
    //     Authorization: `Bearer ${token}`
    //   },
    //   timeout: false
    // });
  },
  async deleteExperiment(groupId: number, experimentId: number) {
    return ApiManager.api.delete(`groups/${groupId}/experiments/${experimentId}`).json<IExperiment>();
  },
  async getExperiment(groupId: number, experimentId: number) {
    return ApiManager.api.get(`groups/${groupId}/experiments/${experimentId}`).json<IExperiment>();
  },
  async getExperimentData(groupId: number, experimentId: number) {
    return ApiManager.api.get(`groups/${groupId}/experiments/${experimentId}/data`).json<IExperimentData>();
  },
  async getChannelStats(acquisitionId: number, channelName: string) {
    const url = `acquisitions/${acquisitionId}/${channelName}/stats?bins=40`;
    let cache;
    if (cacheAvailable) {
      cache = await self.caches.open("stats");
      const found = await cache.match(url);
      if (found) {
        return found.json() as Promise<IChannelStats>;
      }
    }
    const response = await ApiManager.api.get(url);
    if (response.ok) {
      if (cacheAvailable) {
        await cache.put(url, response.clone());
      }
    }
    return response.json() as Promise<IChannelStats>;
  },
  async downloadChannelStackImage(params: IChannelStack) {
    return ApiManager.api.post(`acquisitions/stack`, {
      json: params,
      timeout: false,
    });
  },
  async deleteSlide(id: number) {
    return ApiManager.api.delete(`slides/${id}`).json<ISlide>();
  },
};
