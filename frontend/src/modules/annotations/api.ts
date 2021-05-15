import { ApiManager } from "@/utils/api";

export const api = {
  async deleteGate(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/gates/${id}`).json<number>();
  },
};
