import { ApiManager } from "utils/api";
import { IUserProfile, IUserProfileUpdate } from "./models";

export const api = {
  async getUserProfile() {
    return ApiManager.api.get(`users/profile`).json<IUserProfile>();
  },
  async updateUserProfile(data: IUserProfileUpdate) {
    return ApiManager.api
      .patch(`users/profile`, {
        json: data,
      })
      .json<IUserProfile>();
  },
};
