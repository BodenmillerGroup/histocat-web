import { IUserProfile } from "@/modules/user/models";
import { Mutations } from "vuex-smart-module";
import { MainState } from ".";
import { AppNotification } from "./models";
import { ApiManager } from "@/utils/api";

export class MainMutations extends Mutations<MainState> {
  setToken(payload: string) {
    this.state.token = payload;
    ApiManager.init(payload);
  }

  setLoggedIn(payload: boolean) {
    this.state.isLoggedIn = payload;
  }

  setLogInError(payload: boolean) {
    this.state.logInError = payload;
  }

  setUserProfile(payload: IUserProfile) {
    this.state.userProfile = payload;
  }

  addNotification(payload: AppNotification) {
    this.state.notifications.push(payload);
  }

  removeNotification(payload: AppNotification) {
    this.state.notifications = this.state.notifications.filter((notification) => notification !== payload);
  }

  reset() {
    // acquire initial state
    const s = new MainState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
