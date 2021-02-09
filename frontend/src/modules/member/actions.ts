import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { MemberState } from ".";
import { api } from "./api";
import { MemberGetters } from "./getters";
import { MemberMutations } from "./mutations";
import { IMemberCreate, IMemberUpdate } from "./models";
import { groupModule } from "@/modules/group";

export class MemberActions extends Actions<MemberState, MemberGetters, MemberMutations, MemberActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
  }

  async createMember(payload: IMemberCreate) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.createMember(groupId, payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Group member successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getMember(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getMember(groupId, id);
      if (data) {
        this.mutations.setEntity(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateMember(payload: { id: number; data: IMemberUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateMember(groupId, payload.id, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Group member successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteMember(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteMember(groupId, id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Group member successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getGroupMembers(groupId: number) {
    try {
      const data = await api.getGroupMembers(groupId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
