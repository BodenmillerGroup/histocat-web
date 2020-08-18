import { Module } from "vuex-smart-module";
import { MemberGetters } from "./getters";
import { MemberMutations } from "./mutations";
import { MemberActions } from "./actions";
import { schema } from "normalizr";
import { IMember } from "./models";

export const memberSchema = new schema.Entity("members");
export const memberListSchema = [memberSchema];

export class MemberState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IMember } = {};
}

export const memberModule = new Module({
  namespaced: true,

  state: MemberState,
  getters: MemberGetters,
  mutations: MemberMutations,
  actions: MemberActions,
});
