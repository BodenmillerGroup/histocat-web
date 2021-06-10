import { IDataset, IDatasetUpdate } from "./models";
import { ApiManager } from "@/utils/api";
import { apiUrl } from "@/env";

export const api = {
  async getProjectDatasets(groupId: number, projectId: number) {
    return ApiManager.api.get(`groups/${groupId}/projects/${projectId}/datasets`).json<IDataset[]>();
  },
  async getDataset(groupId: number, datasetId: number) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${datasetId}`).json<IDataset>();
  },
  async updateDataset(groupId: number, datasetId: number, data: IDatasetUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/datasets/${datasetId}`, {
        json: data,
      })
      .json<IDataset>();
  },
  async deleteDataset(groupId: number, datasetId: number) {
    return ApiManager.api.delete(`groups/${groupId}/datasets/${datasetId}`).json();
  },
  async downloadDataset(datasetId: number) {
    return ApiManager.api.get(`datasets/${datasetId}/download`, {
      timeout: false,
    });
  },
  async uploadDataset(
    token: string,
    groupId: number,
    projectId: number,
    formData: FormData,
    onLoadStart: () => void,
    onLoad: () => void,
    onProgress: (event: ProgressEvent) => void,
    onError: () => void
  ) {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", `${apiUrl}/groups/${groupId}/projects/${projectId}/datasets/upload`);
    xhr.setRequestHeader("Authorization", `Bearer ${token}`);

    xhr.upload.onloadstart = onLoadStart;
    xhr.upload.onprogress = onProgress;
    xhr.upload.onload = onLoad;
    xhr.upload.onerror = function () {
      console.log(`Error during file upload: ${xhr.status}.`);
      onError();
    };

    xhr.send(formData);

    // return ky.post(`${apiUrl}/api/v1/projects/${id}/upload`, {
    //   body: data,
    //   headers: {
    //     Authorization: `Bearer ${token}`
    //   },
    //   timeout: false
    // });
  },
};
