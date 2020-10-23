import { apiUrl } from "@/env";
import {
  IChannelStack,
  IChannelStats,
  IChannelUpdate,
  IProject,
  IProjectCreate,
  IProjectData,
  IProjectUpdate,
  ISlide
} from "./models";
import { ApiManager } from "@/utils/api";

const cacheAvailable = false; // 'caches' in self;

export const api = {
  async getGroupProjects(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects`).json<IProject[]>();
  },
  async getProjectsTags(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/tags`).json<string[]>();
  },
  async updateProject(groupId: number, projectId: number, data: IProjectUpdate) {
    return ApiManager.api
      .put(`groups/${groupId}/projects/${projectId}`, {
        json: data,
      })
      .json<IProject>();
  },
  async createProject(data: IProjectCreate) {
    return ApiManager.api
      .post(`groups/${data.group_id}/projects`, {
        json: data,
      })
      .json<IProject>();
  },
  async upload(
    token: string,
    groupId: number,
    projectId: number,
    data,
    onLoadStart: () => void,
    onLoad: () => void,
    onProgress: (event: ProgressEvent) => void,
    onError: () => void
  ) {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", `${apiUrl}/groups/${groupId}/projects/${projectId}/upload`);
    xhr.setRequestHeader("Authorization", `Bearer ${token}`);

    xhr.upload.onloadstart = onLoadStart;
    xhr.upload.onprogress = onProgress;
    xhr.upload.onload = onLoad;
    xhr.upload.onerror = function () {
      console.log(`Error during file upload: ${xhr.status}.`);
      onError();
    };

    xhr.send(data);

    // return ky.post(`${apiUrl}/api/v1/projects/${id}/upload`, {
    //   body: data,
    //   headers: {
    //     Authorization: `Bearer ${token}`
    //   },
    //   timeout: false
    // });
  },
  async deleteProject(groupId: number, projectId: number) {
    return ApiManager.api.delete(`groups/${groupId}/projects/${projectId}`).json<IProject>();
  },
  async getProject(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}`).json<IProject>();
  },
  async getProjectData(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/data`).json<IProjectData>();
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
  async updateChannel(groupId: number, acquisitionId: number, data: IChannelUpdate) {
    return ApiManager.api
      .put(`groups/${groupId}/acquisitions/${acquisitionId}`, {
        json: data,
      })
      .json<IProjectData>();
  },
};
