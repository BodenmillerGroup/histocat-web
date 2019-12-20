import { mainModule } from "@/modules/main";
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
      const data = await api.getUsers();
      if (data) {
        this.mutations.setUsers(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateUser(payload: { id: number; user: IUserProfileUpdate }) {
    try {
      const data = await api.updateUser(payload.id, payload.user);
      this.mutations.setUser(data);
      this.main!.mutations.addNotification({ content: "User successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createUser(payload: IUserProfileCreate) {
    try {
      const data = await api.createUser(payload);
      this.mutations.setUser(data);
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
      const data = await api.signUp(payload);
      this.main!.actions.routeLogOut();
      this.main!.mutations.addNotification({ content: "User successfully signed up", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
