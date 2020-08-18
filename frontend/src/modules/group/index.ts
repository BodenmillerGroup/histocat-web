import { Module } from "vuex-smart-module";
import { GroupGetters } from "./getters";
import { GroupMutations } from "./mutations";
import { GroupActions } from "./actions";
import { schema } from "normalizr";
import { IGroup } from "./models";
import { IMember } from "@/modules/member/models";

export const groupSchema = new schema.Entity("groups");
export const groupListSchema = [groupSchema];

export class GroupState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IGroup } = {};
  activeGroupId: number | null = null;
  myMember: IMember | null = null;
  tags: string[] = [];
}

export const groupModule = new Module({
  namespaced: true,

  state: GroupState,
  getters: GroupGetters,
  mutations: GroupMutations,
  actions: GroupActions,
});
