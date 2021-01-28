import { IModel, IModelUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async updateModel(groupId: number, modelId: number, data: IModelUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/models/${modelId}`, {
        json: data,
      })
      .json<IModel>();
  },
  async getGroupModels(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/models`).json<IModel[]>();
  },
  async getModel(groupId: number, modelId: number) {
    return ApiManager.api.get(`groups/${groupId}/models/${modelId}`).json<IModel>();
  },
  async deleteModel(groupId: number, modelId: number) {
    return ApiManager.api.delete(`groups/${groupId}/models/${modelId}`).json<number>();
  },
  // async createModel(
  //   token: string,
  //   validationId: number,
  //   formData: FormData,
  //   onLoadStart: () => void,
  //   onLoad: () => void,
  //   onProgress: (event: ProgressEvent) => void,
  //   onError: () => void
  // ) {
  //   const xhr = new XMLHttpRequest();
  //   xhr.open("POST", `${apiUrl}/validation/${validationId}/upload`);
  //   xhr.setRequestHeader("Authorization", `Bearer ${token}`);
  //   xhr.upload.onloadstart = onLoadStart;
  //   xhr.upload.onprogress = onProgress;
  //   xhr.upload.onload = onLoad;
  //   xhr.upload.onerror = function() {
  //     console.log(`Error during file upload: ${xhr.status}.`);
  //     onError();
  //   };
  //   xhr.send(formData);
  // },
  async createModel(groupId: number, formData: FormData) {
    return ApiManager.api
      .post(`groups/${groupId}/models`, {
        body: formData,
        timeout: false,
      })
      .json<IModel>();
  },
};
