import ky from "ky";
import { apiUrl } from "@/env";
import { IUserProfile, IUserProfileCreate, IUserProfileUpdate } from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async logInGetToken(username: string, password: string) {
    const params = new URLSearchParams();
    params.append("username", username);
    params.append("password", password);

    return ky.post(`${apiUrl}/auth/login`, { body: params }).json();
  },
  async getMe() {
    return ApiManager.api.get(`users/profile`).json<IUserProfile>();
  },
  async updateMe(data: IUserProfileUpdate) {
    return ApiManager.api
      .patch(`users/profile`, {
        json: data,
      })
      .json<IUserProfile>();
  },
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
  async passwordRecovery(email: string) {
    return ky.post(`${apiUrl}/auth/password-recovery/${email}`).json();
  },
  async checkUserExists(email: string) {
    return ky.get(`${apiUrl}/users/check/${email}`).json<{ exists: boolean }>();
  },
  async signUp(data: IUserProfileCreate) {
    return ky
      .post(`${apiUrl}/auth/signup`, {
        json: data,
      })
      .json();
  },
  async resetPassword(password: string, token: string) {
    return ky
      .post(`${apiUrl}/auth/reset-password/`, {
        json: {
          new_password: password,
          token,
        },
      })
      .json();
  },
};
