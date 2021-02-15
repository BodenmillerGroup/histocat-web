import create from "zustand";
import { IUserProfile, IUserProfileUpdate } from "./models";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "utils/toaster";

type ProfileState = {
  userProfile: IUserProfile | null;

  hasAdminAccess(): boolean | null;
  getUserProfile(): Promise<void>;
  updateUserProfile(params: IUserProfileUpdate): Promise<void>;
};

export const useProfileStore = create<ProfileState>((set, get) => ({
  userProfile: null,

  hasAdminAccess() {
    const userProfile = get().userProfile;
    return userProfile && userProfile.is_admin && userProfile.is_active;
  },

  async getUserProfile() {
    try {
      const data = await api.getUserProfile();
      if (data) {
        set((state) => ({ userProfile: data }));
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateUserProfile(params: IUserProfileUpdate) {
    try {
      const data = await api.updateUserProfile(params);
      set({ userProfile: data });
      AppToaster.show({ message: "Profile successfully updated", intent: "success" });
    } catch (error) {
      displayApiError(error);
    }
  },
}));
