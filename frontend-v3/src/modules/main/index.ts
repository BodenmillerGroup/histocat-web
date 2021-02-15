import create from "zustand";
import { ViewMode } from "./models";

type MainState = {
  showWorkspace: boolean;
  showOptions: boolean;
  viewMode: ViewMode;
  processing: boolean;
  processingProgress: number;

  setProcessing: (payload: boolean) => void;
  setProcessingProgress: (payload: number) => void;
  setViewMode: (payload: ViewMode) => void;
};

export const useMainStore = create<MainState>((set, get) => ({
  showWorkspace: true,
  showOptions: true,
  viewMode: "image",
  processing: false,
  processingProgress: 0,

  setProcessing(payload: boolean) {
    set({ processing: payload });
  },

  setProcessingProgress(payload: number) {
    set({ processingProgress: payload });
  },

  setViewMode(payload: ViewMode) {
    set({ viewMode: payload });
  },
}));
