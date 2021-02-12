import { ApiManager } from "utils/api";
import { IGroup, IGroupCreate, IGroupUpdate } from "./models";
import { IMember } from "../members/models";

export const api = {
  async getGroups() {
    return ApiManager.api.get(`groups`).json<IGroup[]>();
  },
  async getTags() {
    return ApiManager.api.get(`groups/tags`).json<string[]>();
  },
  async createGroup(data: IGroupCreate) {
    return ApiManager.api
      .post(`groups/`, {
        json: data,
      })
      .json<IGroup>();
  },
  async getGroup(id: number) {
    return ApiManager.api.get(`groups/${id}`).json<IGroup>();
  },
  async updateGroup(id: number, data: IGroupUpdate) {
    return ApiManager.api
      .patch(`groups/${id}`, {
        json: data,
      })
      .json<IGroup>();
  },
  async deleteGroup(id: number) {
    return ApiManager.api.delete(`groups/${id}`).json<number>();
  },
  async joinGroup(id: number) {
    return ApiManager.api.post(`groups/${id}/join`).json<IGroup>();
  },
  async getMyMember(groupId: number) {
    return ApiManager.api.get(`groups/${groupId}/members/me`).json<IMember>();
  },
};
