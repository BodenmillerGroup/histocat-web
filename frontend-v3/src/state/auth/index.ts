import create from "zustand";
import { api } from "./api";
import { useNotificationsStore } from "state/notifications";
import { useProfileStore } from "state/profile";
import { ApiManager } from "utils/api";
import { saveLocalToken, getLocalToken, removeLocalToken } from "./helpers";

const addNotification = useNotificationsStore.getState().addNotification;
const getUserProfile = useProfileStore.getState().getUserProfile;

type AuthState = {
  token: string | null;
  loggedIn: boolean | null;
  loginError: boolean;
  login: (username: string, password: string) => void;
  removeLogin: () => void;
  checkLoggedIn: () => void;
};

export const useAuthStore = create<AuthState>((set, get) => ({
  token: "",
  loggedIn: null,
  loginError: false,

  login: async (username: string, password: string) => {
    try {
      const response: any = await api.login(username, password);
      const token = response.access_token;
      if (token) {
        saveLocalToken(token);
        set({ token: token, loggedIn: true, loginError: false });
        ApiManager.init(token);
        await getUserProfile();
        // await this.actions.routeLoggedIn();
        addNotification({ content: "Logged in", color: "success" });
      } else {
        // await this.actions.logOut();
      }
    } catch (error) {
      set({ loggedIn: false, loginError: true });
      // await this.actions.logOut();
    }
  },

  removeLogin: () => {
    removeLocalToken();
    set({ token: "", loggedIn: false });
  },

  checkLoggedIn: async () => {
    if (!get().loggedIn) {
      let token = get().token;
      if (!token) {
        const localToken = getLocalToken();
        if (localToken) {
          set({token: localToken})
          token = localToken;
        }
      }
      if (token) {
        try {
          await getUserProfile();
          set({loggedIn: true})
          ApiManager.init(token);
        } catch (error) {
          get().removeLogin();
        }
      } else {
        get().removeLogin();
      }
    }
  }
}));
