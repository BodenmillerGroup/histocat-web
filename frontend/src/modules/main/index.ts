import { IUserProfile } from "@/modules/user/models";
import { Module } from "vuex-smart-module";
import { MainActions } from "./actions";
import { MainGetters } from "./getters";
import { AppNotification, ViewMode } from "./models";
import { MainMutations } from "./mutations";

export class MainState {
  token = "";
  isLoggedIn: boolean | null = null;
  logInError = false;
  userProfile: IUserProfile | null = null;
  dashboardMiniDrawer = true;
  dashboardShowDrawer = true;
  notifications: AppNotification[] = [];
  showWorkspace = true;
  showOptions = true;
  viewMode: ViewMode = "image";

  processing = false;
  processingProgress = 0;

  maskMode: "raw" | "mask" | "origin" = "raw";
  mouseMode: "panZoom" | "lasso" | "rotate" = "panZoom";
}

export const mainModule = new Module({
  namespaced: true,

  state: MainState,
  getters: MainGetters,
  mutations: MainMutations,
  actions: MainActions,
});
