import ky from "ky";
import { apiUrl } from "env";

export const api = {
  async login(username: string, password: string) {
    const params = new URLSearchParams();
    params.append("username", username);
    params.append("password", password);

    return ky.post(`${apiUrl}/auth/login`, { body: params }).json();
  },
};
