import { Getters } from "vuex-smart-module";
import { UiState } from ".";

export class UiGetters extends Getters<UiState> {
  get activeLayout() {
    return this.state.activeLayout;
  }

  get responsive() {
    return this.state.responsive;
  }

  get dashboardShowDrawer() {
    return this.state.dashboardShowDrawer;
  }

  get dashboardMiniDrawer() {
    return this.state.dashboardMiniDrawer;
  }

  get showWorkspace() {
    return this.state.showWorkspace;
  }

  get showOptions() {
    return this.state.showOptions;
  }

  get processing() {
    return this.state.processing;
  }

  get processingProgress() {
    return this.state.processingProgress;
  }

  get viewMode() {
    return this.state.viewMode;
  }

  get maskMode() {
    return this.state.maskMode;
  }

  get maskOpacity() {
    return this.state.maskOpacity;
  }

  get mouseMode() {
    return this.state.mouseMode;
  }
}
