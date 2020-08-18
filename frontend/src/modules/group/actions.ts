import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GroupState } from ".";
import { api } from "./api";
import { GroupGetters } from "./getters";
import { GroupMutations } from "./mutations";
import { saveAs } from "file-saver";
import { IGroupCreate, IGroupUpdate } from "./models";

export class GroupActions extends Actions<GroupState, GroupGetters, GroupMutations, GroupActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
  }

  async getGroups() {
    try {
      const data = await api.getGroups();
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createGroup(payload: IGroupCreate) {
    try {
      const data = await api.createGroup(payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Group successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getGroup(id: number) {
    try {
      const data = await api.getGroup(id);
      if (data) {
        this.mutations.setEntity(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateGroup(payload: { id: number; data: IGroupUpdate }) {
    try {
      const data = await api.updateGroup(payload.id, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Group successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteGroup(id: number) {
    try {
      const data = await api.deleteGroup(id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Group successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async exportGroupData(payload: { id: number; format: "json" | "csv" }) {
    try {
      const blob = await api.exportGroupData(payload.id, payload.format);
      saveAs(blob, `group_${payload.id}.zip`);
      this.main!.mutations.addNotification({ content: "Group data successfully exported", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async importGroupData(payload: { formData: FormData }) {
    try {
      const data = await api.importGroupData(payload.formData);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Group data successfully imported", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async joinGroup(id: number) {
    try {
      const result = await api.joinGroup(id);
      if (result) {
        this.main!.mutations.addNotification({
          content: "Request to join the group successfully submitted",
          color: "success",
        });
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getMyMember(groupId: number) {
    try {
      const data = await api.getMyMember(groupId);
      if (data) {
        this.mutations.setMyMember(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
