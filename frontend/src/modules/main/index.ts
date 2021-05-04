import { IUserProfile } from "@/modules/user/models";
import { Module } from "vuex-smart-module";
import { MainActions } from "./actions";
import { MainGetters } from "./getters";
import { AppNotification } from "./models";
import { MainMutations } from "./mutations";

export class MainState {
  token = "";
  isLoggedIn: boolean | null = null;
  logInError = false;
  userProfile: IUserProfile | null = null;
  notifications: AppNotification[] = [];
}

export const mainModule = new Module({
  namespaced: true,

  state: MainState,
  getters: MainGetters,
  mutations: MainMutations,
  actions: MainActions,
});
