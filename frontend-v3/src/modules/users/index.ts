import create from "zustand";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "../../utils/toaster";
import { IUserProfile, IUserProfileCreate, IUserProfileUpdate } from "../profile/models";

type UsersState = {
  users: IUserProfile[];

  getUsers(): Promise<void>;
  updateUser(id: number, params: IUserProfileUpdate): Promise<void>;
  createUser(params: IUserProfileCreate): Promise<void>;
};

export const useUsersStore = create<UsersState>((set, get) => ({
  users: [],

  async getUsers() {
    try {
      const data = await api.getUsers();
      if (data) {
        set({ users: data });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateUser(id: number, params: IUserProfileUpdate) {
    try {
      const data = await api.updateUser(id, params);
      if (data) {
        const users = get().users.filter((user: IUserProfile) => user.id !== data.id);
        users.push(data);
        set({ users: users });
        AppToaster.show({ message: "User successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async createUser(params: IUserProfileCreate) {
    try {
      const data = await api.createUser(params);
      if (data) {
        const users = get().users.filter((user: IUserProfile) => user.id !== data.id);
        users.push(data);
        set({ users: users });
        AppToaster.show({ message: "User successfully created", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },
}));
