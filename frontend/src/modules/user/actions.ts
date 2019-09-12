import { mainModule } from "@/modules/main";
import router from "@/router";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { UserState } from ".";
import { api } from "./api";
import { UserGetters } from "./getters";
import { IUserProfileCreate, IUserProfileUpdate } from "./models";
import { UserMutations } from "./mutations";

export class UserActions extends Actions<UserState, UserGetters, UserMutations, UserActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
  }

  async getUsers() {
    try {
      const data = await api.getUsers(this.main!.getters.token);
      if (data) {
        this.mutations.setUsers(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateUser(payload: { id: number; user: IUserProfileUpdate }) {
    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.updateUser(this.main!.getters.token, payload.id, payload.user),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setUser(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "User successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createUser(payload: IUserProfileCreate) {
    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createUser(this.main!.getters.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setUser(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "User successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async checkUserExists(email: string) {
    try {
      const data = await api.checkUserExists(email);
      return data.exists;
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async signUp(payload: IUserProfileCreate) {
    try {
      const notification = { content: "signing up", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.signUp(payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.main!.actions.routeLogOut();
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "User successfully signed up", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
