import ky from "ky";
import { apiUrl } from "env";
import { IUserProfileCreate } from "../profile/models";

export const api = {
  async login(username: string, password: string) {
    const params = new URLSearchParams();
    params.append("username", username);
    params.append("password", password);

    return ky.post(`${apiUrl}/auth/login`, { body: params }).json();
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
