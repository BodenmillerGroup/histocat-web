import { IModel, IModelUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async updateModel(modelId: number, data: IModelUpdate) {
    return ApiManager.api
      .patch(`models/${modelId}`, {
        json: data,
      })
      .json<IModel>();
  },
  async getAllModels() {
    return ApiManager.api.get(`models`).json<IModel[]>();
  },
  async getModel(modelId: number) {
    return ApiManager.api.get(`models/${modelId}`).json<IModel>();
  },
  async deleteModel(modelId: number) {
    return ApiManager.api.delete(`models/${modelId}`).json<number>();
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
  async createModel(formData: FormData) {
    return ApiManager.api
      .post(`models`, {
        body: formData,
        timeout: false,
      })
      .json<IModel>();
  },
};
