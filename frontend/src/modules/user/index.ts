import { UserActions } from "@/modules/user/actions";
import { UserMutations } from "@/modules/user/mutations";
import { Module } from "vuex-smart-module";
import { UserGetters } from "./getters";
import { IUserProfile } from "./models";

export class UserState {
  users: IUserProfile[] = [];
}

export const userModule = new Module({
  namespaced: false,

  state: UserState,
  getters: UserGetters,
  mutations: UserMutations,
  actions: UserActions
});
