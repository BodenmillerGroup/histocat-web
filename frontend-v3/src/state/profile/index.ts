import create from "zustand";
import { IUserProfile } from "./models";
import { api } from "./api";
import { useNotificationsStore } from "../notifications";

const displayApiError = useNotificationsStore.getState().displayApiError;

type ProfileState = {
  userProfile: IUserProfile | null;
  getUserProfile: () => void;
};

export const useProfileStore = create<ProfileState>((set, get) => ({
  userProfile: null,

  getUserProfile: async () => {
    try {
      const data = await api.getUserProfile();
      if (data) {
        set((state) => ({ userProfile: data }));
      }
    } catch (error) {
      displayApiError(error);
    }
  },
}));
