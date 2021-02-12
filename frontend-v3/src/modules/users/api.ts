import { ApiManager } from "utils/api";
import { IUserProfile, IUserProfileCreate, IUserProfileUpdate } from "../profile/models";

export const api = {
  async getUsers() {
    return ApiManager.api.get(`users`).json<IUserProfile[]>();
  },
  async updateUser(id: number, data: IUserProfileUpdate) {
    return ApiManager.api
      .put(`users/${id}`, {
        json: data,
      })
      .json<IUserProfile>();
  },
  async createUser(data: IUserProfileCreate) {
    return ApiManager.api
      .post(`users`, {
        json: data,
      })
      .json<IUserProfile>();
  },
};
