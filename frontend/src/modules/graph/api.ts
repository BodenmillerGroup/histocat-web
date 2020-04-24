import { ApiManager } from "@/utils/api";

export const api = {
  async getSchema() {
    return ApiManager.api.get(`schema`).json<any>();
  },

  async getObsAnnotation(name: string) {
    const url = `annotations/obs?annotation_name=${encodeURIComponent(name)}`;
    return ApiManager.api
      .get(url, {
        headers: {
          "content-type": "application/octet-stream",
        },
      })
      .arrayBuffer();
  },

  async getVarAnnotation(name: string) {
    const url = `annotations/var?annotation_name=${encodeURIComponent(name)}`;
    return ApiManager.api
      .get(url, {
        headers: {
          "content-type": "application/octet-stream",
        },
      })
      .arrayBuffer();
  },
};
