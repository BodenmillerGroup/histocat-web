import { settingsModule } from "@/modules/settings";
import { api } from "@/modules/user/api";
import router from "@/router";
import { getLocalToken, removeLocalToken, saveLocalToken } from "@/utils/auth";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { MainState } from ".";
import { MainGetters } from "./getters";
import { AppNotification } from "./models";
import { MainMutations } from "./mutations";

export class MainActions extends Actions<MainState, MainGetters, MainMutations, MainActions> {
  // Declare context type
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.settings = settingsModule.context(store);
  }

  async logIn(payload: { username: string; password: string }) {
    try {
      const data: any = await api.logInGetToken(payload.username, payload.password);
      const token = data.access_token;
      if (token) {
        saveLocalToken(token);
        this.mutations.setToken(token);
        this.mutations.setLoggedIn(true);
        this.mutations.setLogInError(false);
        await this.actions.getUserProfile();
        await this.actions.routeLoggedIn();
        this.mutations.addNotification({ content: "Logged in", color: "success" });
      } else {
        await this.actions.logOut();
      }
    } catch (err) {
      this.mutations.setLogInError(true);
      await this.actions.logOut();
    }
  }

  async getUserProfile() {
    try {
      const data = await api.getMe();
      if (data) {
        this.mutations.setUserProfile(data);
      }
    } catch (error) {
      await this.actions.checkApiError(error);
    }
  }

  async updateUserProfile(payload) {
    try {
      const data = await api.updateMe(payload);
      this.mutations.setUserProfile(data);
      this.mutations.addNotification({ content: "Profile successfully updated", color: "success" });
    } catch (error) {
      await this.actions.checkApiError(error);
    }
  }

  async checkLoggedIn() {
    if (!this.getters.isLoggedIn) {
      let token = this.getters.token;
      if (!token) {
        const localToken = getLocalToken();
        if (localToken) {
          this.mutations.setToken(localToken);
          token = localToken;
        }
      }
      if (token) {
        try {
          const data = await api.getMe();
          this.mutations.setLoggedIn(true);
          this.mutations.setUserProfile(data);
        } catch (error) {
          await this.actions.removeLogIn();
        }
      } else {
        await this.actions.removeLogIn();
      }
    }
  }

  async removeLogIn() {
    removeLocalToken();
    this.mutations.setToken("");
    this.mutations.setLoggedIn(false);
  }

  async logOut() {
    await this.actions.removeLogIn();
    await this.actions.routeLogOut();
  }

  async userLogOut() {
    await this.actions.logOut();
    this.mutations.addNotification({ content: "Logged out", color: "success" });
  }

  routeLogOut() {
    if (router.currentRoute.path !== "/login") {
      router.push("/login", () => {});
    }
  }

  async checkApiError(error) {
    this.mutations.addNotification({ content: error.message, color: "error" });
    if (error.response) {
      console.error("API error: ", error.response);
      if (error.response.status === 401) {
        // await this.actions.logOut();
      }
    }
  }

  routeLoggedIn() {
    if (router.currentRoute.path === "/login" || router.currentRoute.path === "/") {
      router.push("/main", () => {});
    }
  }

  async removeNotification(payload: { notification: AppNotification; timeout: number }) {
    return new Promise((resolve, reject) => {
      setTimeout(() => {
        this.mutations.removeNotification(payload.notification);
        resolve(true);
      }, payload.timeout);
    });
  }

  async passwordRecovery(payload: { username: string }) {
    try {
      await api.passwordRecovery(payload.username);
      this.mutations.addNotification({ content: "Password recovery email sent", color: "success" });
      await this.actions.logOut();
    } catch (error) {
      await this.actions.checkApiError(error);
    }
  }

  async resetPassword(payload: { password: string; token: string }) {
    try {
      const response = await api.resetPassword(payload.password, payload.token);
      this.mutations.addNotification({ content: "Password successfully reset", color: "success" });
      await this.actions.logOut();
    } catch (error) {
      await this.actions.checkApiError(error);
    }
  }
}
