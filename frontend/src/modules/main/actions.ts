import { settingsModule } from "@/modules/settings";
import { api } from "@/modules/user/api";
import router from "@/router";
import { getLocalToken, removeLocalToken, saveLocalToken } from "@/utils";
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
        await this.getUserProfile();
        await this.routeLoggedIn();
        this.mutations.addNotification({ content: "Logged in", color: "success" });
      } else {
        await this.logOut();
      }
    } catch (err) {
      this.mutations.setLogInError(true);
      await this.logOut();
    }
  }

  async getUserProfile() {
    try {
      const data = await api.getMe(this.state.token);
      if (data) {
        this.mutations.setUserProfile(data);
      }
    } catch (error) {
      await this.checkApiError(error);
    }
  }

  async updateUserProfile(payload) {
    try {
      const loadingNotification = { content: "saving", showProgress: true };
      this.mutations.addNotification(loadingNotification);
      const data = (await Promise.all([
        api.updateMe(this.state.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setUserProfile(data);
      this.mutations.removeNotification(loadingNotification);
      this.mutations.addNotification({ content: "Profile successfully updated", color: "success" });
    } catch (error) {
      await this.checkApiError(error);
    }
  }

  async checkLoggedIn() {
    if (!this.state.isLoggedIn) {
      let token = this.state.token;
      if (!token) {
        const localToken = getLocalToken();
        if (localToken) {
          this.mutations.setToken(localToken);
          token = localToken;
        }
      }
      if (token) {
        try {
          const data = await api.getMe(token);
          this.mutations.setLoggedIn(true);
          this.mutations.setUserProfile(data);
        } catch (error) {
          await this.removeLogIn();
        }
      } else {
        await this.removeLogIn();
      }
    }
  }

  async removeLogIn() {
    removeLocalToken();
    this.mutations.setToken("");
    this.mutations.setLoggedIn(false);
  }

  async logOut() {
    await this.removeLogIn();
    await this.routeLogOut();
  }

  async userLogOut() {
    await this.logOut();
    this.mutations.addNotification({ content: "Logged out", color: "success" });
  }

  routeLogOut() {
    if (router.currentRoute.path !== "/login") {
      router.push("/login");
    }
  }

  async checkApiError(payload) {
    console.log("API error: ", payload);
    if (payload.response && payload.response.status === 401) {
      await this.logOut();
    }
  }

  routeLoggedIn() {
    if (router.currentRoute.path === "/login" || router.currentRoute.path === "/") {
      router.push("/main");
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
    const loadingNotification = { content: "Sending password recovery email", showProgress: true };
    try {
      this.mutations.addNotification(loadingNotification);
      const response = (await Promise.all([
        api.passwordRecovery(payload.username),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.removeNotification(loadingNotification);
      this.mutations.addNotification({ content: "Password recovery email sent", color: "success" });
      await this.logOut();
    } catch (error) {
      this.mutations.removeNotification(loadingNotification);
      this.mutations.addNotification({ color: "error", content: "Incorrect username" });
    }
  }

  async resetPassword(payload: { password: string; token: string }) {
    const loadingNotification = { content: "Resetting password", showProgress: true };
    try {
      this.mutations.addNotification(loadingNotification);
      const response = (await Promise.all([
        api.resetPassword(payload.password, payload.token),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.removeNotification(loadingNotification);
      this.mutations.addNotification({ content: "Password successfully reset", color: "success" });
      await this.logOut();
    } catch (error) {
      this.mutations.removeNotification(loadingNotification);
      this.mutations.addNotification({ color: "error", content: "Error resetting password" });
    }
  }
}
