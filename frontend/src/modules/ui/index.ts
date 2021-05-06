import { Module } from "vuex-smart-module";
import { UiActions } from "./actions";
import { UiGetters } from "./getters";
import { UiMutations } from "./mutations";
import { GoldenLayout } from "golden-layout";
import { IResponsive } from "./models";
import { DEFAULT_LAYOUTS } from "./defaultLayouts";

export const PROJECT_LAYOUTS_STORAGE_KEY = "ProjectLayouts";

let initialLayouts = DEFAULT_LAYOUTS;
if (localStorage.getItem(PROJECT_LAYOUTS_STORAGE_KEY)) {
  initialLayouts = JSON.parse(localStorage.getItem(PROJECT_LAYOUTS_STORAGE_KEY)!);
}

export class UiState {
  layouts = initialLayouts;
  activeLayout = initialLayouts[0];
  goldenLayout: GoldenLayout | null = null;

  responsive: IResponsive = {
    width: null,
    height: null,
  };

  processing = false;
  processingProgress = 0;

  maskMode: "raw" | "mask" | "origin" = "raw";
  maskOpacity = 1.0;
  mouseMode: "panZoom" | "lasso" | "rotate" = "panZoom";
}

export const uiModule = new Module({
  namespaced: true,

  state: UiState,
  getters: UiGetters,
  mutations: UiMutations,
  actions: UiActions,
});
