import { apiUrl } from "@/env";
import {
  IChannelStack,
  IChannelStats,
  IChannelUpdate,
  IProject,
  IProjectCreate,
  IProjectData,
  IProjectUpdate,
  ISlide,
} from "./models";
import { ApiManager } from "@/utils/api";

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
  async uploadSlide(
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
    xhr.open("POST", `${apiUrl}/groups/${groupId}/projects/${projectId}/slides/upload`);
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
  async getChannelStats(groupId: number, acquisitionId: number, channelName: string) {
    const url = `groups/${groupId}/acquisitions/${acquisitionId}/${channelName}/stats?bins=40`;
    const response = await ApiManager.api.get(url);
    return response.json() as Promise<IChannelStats>;
  },
  async downloadChannelStackImage(groupId: number, params: IChannelStack) {
    return ApiManager.api.post(`groups/${groupId}/acquisitions/stack`, {
      json: params,
      timeout: false,
    });
  },
  async downloadOmeTiffImage(groupId: number, acquisitionId: number) {
    return ApiManager.api.get(`groups/${groupId}/acquisitions/${acquisitionId}/download`, {
      timeout: false,
    });
  },
  async deleteSlide(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/slides/${id}`).json<ISlide>();
  },
  async updateChannel(groupId: number, acquisitionId: number, data: IChannelUpdate) {
    return ApiManager.api
      .put(`groups/${groupId}/acquisitions/${acquisitionId}`, {
        json: data,
      })
      .json<IProjectData>();
  },
};
