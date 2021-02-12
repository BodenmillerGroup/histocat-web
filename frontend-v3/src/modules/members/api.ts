import { ApiManager } from "utils/api";
import { IMember, IMemberCreate, IMemberUpdate } from "./models";

export const api = {
  async createMember(groupId: number, data: IMemberCreate) {
    return ApiManager.api
      .post(`groups/${groupId}/members`, {
        json: data,
      })
      .json<IMember>();
  },
  async getMember(groupId: number, id: number) {
    return ApiManager.api.get(`groups/${groupId}/members/${id}`).json<IMember>();
  },
  async updateMember(groupId: number, id: number, data: IMemberUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/members/${id}`, {
        json: data,
      })
      .json<IMember>();
  },
  async deleteMember(groupId: number, id: number) {
    return ApiManager.api.delete(`groups/${groupId}/members/${id}`).json<number>();
  },
  async getGroupMembers(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/members`).json<IMember[]>();
  },
};
