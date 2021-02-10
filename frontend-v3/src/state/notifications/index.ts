import create from "zustand";
import { AppNotification } from "./models";

type NotificationsState = {
  notifications: AppNotification[];
  addNotification: (notification: AppNotification) => void;
  removeNotification: (notification: AppNotification) => void;
  displayApiError: (error: any) => void;
};

export const useNotificationsStore = create<NotificationsState>((set, get) => ({
  notifications: [],

  addNotification: (notification: AppNotification) => {
    set((state) => ({ notifications: state.notifications.concat(notification) }));
  },

  removeNotification: (notification: AppNotification) => {
    set((state) => ({ notifications: state.notifications.filter((item) => item !== notification) }));
  },

  displayApiError: (error: any) => {
    get().addNotification({ content: error.message, color: "error" });
    if (error.response) {
      console.error("API error: ", error.response);
      if (error.response.status === 401) {
        // await this.actions.logOut();
      }
    }
  },
}));
