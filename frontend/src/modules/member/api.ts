import { ApiManager } from "@/utils/api";
import { IMember, IMemberCreate, IMemberUpdate } from "./models";

export const api = {
  async createMember(data: IMemberCreate) {
    return ApiManager.api
      .post(`members`, {
        json: data,
      })
      .json<IMember>();
  },
  async getMember(id: number) {
    return ApiManager.api.get(`members/${id}`).json<IMember>();
  },
  async updateMember(id: number, data: IMemberUpdate) {
    return ApiManager.api
      .patch(`members/${id}`, {
        json: data,
      })
      .json<IMember>();
  },
  async deleteMember(id: number) {
    return ApiManager.api.delete(`members/${id}`).json<number>();
  },
  async getGroupMembers(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/members`).json<IMember[]>();
  },
};
