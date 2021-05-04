import { Mutations } from "vuex-smart-module";
import { UiState } from ".";
import { ILayout, IResponsive, ViewMode } from "./models";

export class UiMutations extends Mutations<UiState> {
  setActiveLayout(value: ILayout) {
    this.state.activeLayout = value;
  }

  setResponsive(value: IResponsive) {
    this.state.responsive = value;
  }

  setDashboardMiniDrawer(payload: boolean) {
    this.state.dashboardMiniDrawer = payload;
  }

  setDashboardShowDrawer(payload: boolean) {
    this.state.dashboardShowDrawer = payload;
  }

  setLayout(payload: { showWorkspace: boolean; showOptions: boolean }) {
    this.state.showWorkspace = payload.showWorkspace;
    this.state.showOptions = payload.showOptions;
  }

  setProcessing(payload: boolean) {
    this.state.processing = payload;
  }

  setProcessingProgress(payload: number) {
    this.state.processingProgress = payload;
  }

  setViewMode(value: ViewMode) {
    this.state.viewMode = value;
  }

  setMaskMode(payload: "raw" | "mask" | "origin") {
    this.state.maskMode = payload;
  }

  setMaskOpacity(payload: number) {
    this.state.maskOpacity = payload;
  }

  setMouseMode(mode: "panZoom" | "lasso" | "rotate") {
    this.state.mouseMode = mode;
  }

  reset() {
    // acquire initial state
    const s = new UiState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
