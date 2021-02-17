import create from "zustand";
import { api } from "./api";
import { useProfileStore } from "modules/profile";
import { ApiManager, displayApiError } from "utils/api";
import { saveLocalToken, getLocalToken, removeLocalToken } from "./helpers";
import history from "utils/history";
import { AppToaster } from "utils/toaster";
import { IUserProfileCreate } from "../profile/models";

const getUserProfile = useProfileStore.getState().getUserProfile;

type AuthState = {
  token: string | null;
  loggedIn: boolean | null;
  loginError: boolean;

  login: (username: string, password: string) => Promise<void>;
  checkLoggedIn: () => Promise<void>;
  removeLogin: () => void;
  routeLoggedIn: () => void;
  routeLogout: () => void;
  logout: () => void;
  userLogout: () => void;
  checkUserExists: (email: string) => Promise<boolean>;
  signUp: (payload: IUserProfileCreate) => Promise<void>;
  passwordRecovery: (email: string) => Promise<void>;
  resetPassword: (password: string, token: string) => Promise<void>;
  updateUserPassword: (password: string) => Promise<void>;
};

export const useAuthStore = create<AuthState>((set, get) => ({
  token: null,
  loggedIn: null,
  loginError: false,

  async login(username: string, password: string) {
    try {
      const response: any = await api.login(username, password);
      const token = response.access_token;
      if (token) {
        saveLocalToken(token);
        set({ token: token, loggedIn: true, loginError: false });
        ApiManager.init(token);
        await getUserProfile();
        get().routeLoggedIn();
        AppToaster.show({ message: "Logged in", intent: "success" });
      } else {
        get().logout();
      }
    } catch (error) {
      set({ loginError: true });
      get().logout();
    }
  },

  async checkLoggedIn() {
    if (!get().loggedIn) {
      let token = get().token;
      if (!token) {
        const localToken = getLocalToken();
        if (localToken) {
          set({ token: localToken });
          token = localToken;
        }
      }
      if (token) {
        try {
          ApiManager.init(token);
          set({ loggedIn: true });
          await getUserProfile();
        } catch (error) {
          get().removeLogin();
        }
      } else {
        get().removeLogin();
      }
    }
  },

  removeLogin() {
    removeLocalToken();
    set({ token: null, loggedIn: false });
  },

  routeLoggedIn() {
    if (history.location.pathname === "/login" || history.location.pathname === "/") {
      history.push("/");
    }
  },

  routeLogout() {
    if (history.location.pathname !== "/login") {
      history.push("/login");
    }
  },

  logout() {
    get().routeLogout();
    get().removeLogin();
  },

  userLogout() {
    get().logout();
    AppToaster.show({ message: "Logged out", intent: "success" });
  },

  async checkUserExists(email: string) {
    try {
      const data = await api.checkUserExists(email);
      return data.exists;
    } catch (error) {
      displayApiError(error);
      return false;
    }
  },

  async signUp(params: IUserProfileCreate) {
    try {
      const data = await api.signUp(params);
      get().routeLogout();
      AppToaster.show({ message: "Registration confirmation email was sent", intent: "success" });
    } catch (error) {
      displayApiError(error);
    }
  },

  async passwordRecovery(email: string) {
    try {
      await api.passwordRecovery(email);
      AppToaster.show({ message: "Password recovery email sent", intent: "success" });
      get().logout();
    } catch (error) {
      displayApiError(error);
    }
  },

  async resetPassword(password: string, token: string) {
    try {
      const response = await api.resetPassword(password, token);
      AppToaster.show({ message: "Password successfully reset", intent: "success" });
      get().logout();
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateUserPassword(password: string) {
    try {
      // const data = await api.updateUserProfile(payload);
      get().logout();
      AppToaster.show({ message: "Not implemented", intent: "warning" });
    } catch (error) {
      displayApiError(error);
    }
  },
}));
