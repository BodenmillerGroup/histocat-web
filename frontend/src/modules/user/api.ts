import ky from "ky";
import { apiUrl } from "@/env";
import { IUserProfile, IUserProfileCreate, IUserProfileUpdate } from "./models";

export const api = {
  async logInGetToken(username: string, password: string) {
    const params = new URLSearchParams();
    params.append("username", username);
    params.append("password", password);

    return ky.post(`${apiUrl}/api/v1/auth/access-token`, { body: params }).json();
  },
  async getMe(token: string) {
    return ky
      .get(`${apiUrl}/api/v1/users/me`, {
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IUserProfile>();
  },
  async updateMe(token: string, data: IUserProfileUpdate) {
    return ky
      .put(`${apiUrl}/api/v1/users/me`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IUserProfile>();
  },
  async getUsers(token: string) {
    return ky
      .get(`${apiUrl}/api/v1/users/`, {
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IUserProfile[]>();
  },
  async updateUser(token: string, id: number, data: IUserProfileUpdate) {
    return ky
      .put(`${apiUrl}/api/v1/users/${id}`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IUserProfile>();
  },
  async createUser(token: string, data: IUserProfileCreate) {
    return ky
      .post(`${apiUrl}/api/v1/users/`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IUserProfile>();
  },
  async passwordRecovery(email: string) {
    return ky.post(`${apiUrl}/api/v1/auth/password-recovery/${email}`).json();
  },
  async checkUserExists(email: string) {
    return ky.get(`${apiUrl}/api/v1/users/check/${email}`).json<{ exists: boolean }>();
  },
  async signUp(data: IUserProfileCreate) {
    return ky
      .post(`${apiUrl}/api/v1/users/signup`, {
        json: data
      })
      .json();
  },
  async resetPassword(password: string, token: string) {
    return ky
      .post(`${apiUrl}/api/v1/auth/reset-password/`, {
        json: {
          new_password: password,
          token
        }
      })
      .json();
  }
};
