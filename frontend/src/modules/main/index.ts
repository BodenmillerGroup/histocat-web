import { IUserProfile } from "@/modules/user/models";
import { Module } from "vuex-smart-module";
import { MainActions } from "./actions";
import { MainGetters } from "./getters";
import { AppNotification } from "./models";
import { MainMutations } from "./mutations";

export class MainState {
  token: string = "";
  isLoggedIn: boolean | null = null;
  logInError: boolean = false;
  userProfile: IUserProfile | null = null;
  dashboardMiniDrawer: boolean = true;
  dashboardShowDrawer: boolean = true;
  notifications: AppNotification[] = [];
  showWorkspace: boolean = true;
  showOptions: boolean = true;

  processing: boolean = false;
  processingProgress: number = 0;
}

export const mainModule = new Module({
  namespaced: false,

  state: MainState,
  getters: MainGetters,
  mutations: MainMutations,
  actions: MainActions
});
