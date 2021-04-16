import { IUserProfile } from "@/modules/user/models";
import { Mutations } from "vuex-smart-module";
import { MainState } from ".";
import { AppNotification, ViewMode } from "./models";
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

  setDashboardMiniDrawer(payload: boolean) {
    this.state.dashboardMiniDrawer = payload;
  }

  setDashboardShowDrawer(payload: boolean) {
    this.state.dashboardShowDrawer = payload;
  }

  addNotification(payload: AppNotification) {
    this.state.notifications.push(payload);
  }

  removeNotification(payload: AppNotification) {
    this.state.notifications = this.state.notifications.filter((notification) => notification !== payload);
  }

  setLayout(payload: { showWorkspace: boolean; showOptions: boolean }) {
    this.state.showWorkspace = payload.showWorkspace;
    this.state.showOptions = payload.showOptions;
  }

  setProcessing(payload: boolean) {
    this.state.processing = payload;
  }

  setProcessingProgress(payload: number) {
    this.state.processingProgress = payload;
  }

  setViewMode(value: ViewMode) {
    this.state.viewMode = value;
  }

  setMaskMode(payload: "raw" | "mask" | "origin") {
    this.state.maskMode = payload;
  }

  setMouseMode(mode: "panZoom" | "lasso" | "rotate") {
    this.state.mouseMode = mode;
  }

  reset() {
    // acquire initial state
    const s = new MainState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
