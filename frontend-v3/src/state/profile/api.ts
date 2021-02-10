import { ApiManager } from "utils/api";
import { IUserProfile } from "./models";

export const api = {
  async getUserProfile() {
    return ApiManager.api.get(`users/profile`).json<IUserProfile>();
  },
};
