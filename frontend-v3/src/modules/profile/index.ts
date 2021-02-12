import create from "zustand";
import { IUserProfile, IUserProfileUpdate } from "./models";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "utils/toaster";

type ProfileState = {
  userProfile: IUserProfile | null;
  getUserProfile: () => Promise<void>;
  updateUserProfile: (payload: IUserProfileUpdate) => Promise<void>;
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

  updateUserProfile: async (payload: IUserProfileUpdate) => {
    try {
      const data = await api.updateUserProfile(payload);
      set({ userProfile: data });
      AppToaster.show({ message: "Profile successfully updated", intent: "success" });
    } catch (error) {
      displayApiError(error);
    }
  },
}));
