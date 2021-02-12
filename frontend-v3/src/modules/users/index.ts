import create from "zustand";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "../../utils/toaster";
import { IUserProfile, IUserProfileCreate, IUserProfileUpdate } from "../profile/models";

type UsersState = {
  users: IUserProfile[];
};

export const useUsersStore = create<UsersState>((set, get) => ({
  users: [],

  getUsers: async () => {
    try {
      const data = await api.getUsers();
      if (data) {
        set({ users: data });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  updateUser: async (payload: { id: number; user: IUserProfileUpdate }) => {
    try {
      const data = await api.updateUser(payload.id, payload.user);
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

  createUser: async (payload: IUserProfileCreate) => {
    try {
      const data = await api.createUser(payload);
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
