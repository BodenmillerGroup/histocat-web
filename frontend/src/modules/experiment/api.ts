import { apiUrl } from "@/env";
import {
  IChannelStack,
  IChannelStats,
  IExperiment,
  IExperimentCreate,
  IExperimentUpdate,
  IShare,
  IShareCreate,
  ISlide
} from "./models";
import { ApiManager } from "@/utils/api";

const cacheAvailable = false; // 'caches' in self;

export const api = {
  async getExperiments() {
    return ApiManager.api.get(`experiments`).json<IExperiment[]>();
  },
  async getTags() {
    return ApiManager.api.get(`experiments/tags`).json<string[]>();
  },
  async updateExperiment(id: number, data: IExperimentUpdate) {
    return ApiManager.api
      .put(`experiments/${id}`, {
        json: data
      })
      .json<IExperiment>();
  },
  async createExperiment(data: IExperimentCreate) {
    return ApiManager.api
      .post(`experiments`, {
        json: data
      })
      .json<IExperiment>();
  },
  async upload(
    token: string,
    id: number,
    data,
    onLoadStart: () => void,
    onLoad: () => void,
    onProgress: (event: ProgressEvent) => void,
    onError: () => void
  ) {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", `${apiUrl}/experiments/${id}/upload`);
    xhr.setRequestHeader("Authorization", `Bearer ${token}`);

    xhr.upload.onloadstart = onLoadStart;
    xhr.upload.onprogress = onProgress;
    xhr.upload.onload = onLoad;
    xhr.upload.onerror = function() {
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
  async deleteExperiment(id: number) {
    return ApiManager.api.delete(`experiments/${id}`).json<IExperiment>();
  },
  async getExperiment(id: number) {
    return ApiManager.api.get(`experiments/${id}`).json<IExperiment>();
  },
  async getExperimentData(id: number) {
    return ApiManager.api.get(`experiments/${id}/data`).json<IExperiment>();
  },
  async getChannelStats(id: number) {
    const url = `channels/${id}/stats?bins=100`;
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
    return ApiManager.api.post(`channels/stack`, {
      json: params,
      timeout: false
    });
  },
  async createShare(data: IShareCreate) {
    return ApiManager.api
      .post(`share`, {
        json: data
      })
      .json<IShare>();
  },
  async getExperimentShares(id: number) {
    return ApiManager.api.get(`share/${id}`).json<IShare[]>();
  },
  async deleteSlide(id: number) {
    return ApiManager.api.delete(`slides/${id}`).json<ISlide>();
  }
};
