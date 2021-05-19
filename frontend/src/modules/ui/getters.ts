import { Getters } from "vuex-smart-module";
import { UiState } from ".";

export class UiGetters extends Getters<UiState> {
  get layouts() {
    return this.state.layouts;
  }

  get activeLayout() {
    return this.state.activeLayout;
  }

  get responsive() {
    return this.state.responsive;
  }

  get processing() {
    return this.state.processing;
  }

  get processingProgress() {
    return this.state.processingProgress;
  }

  get showMask() {
    return this.state.showMask;
  }

  get maskOpacity() {
    return this.state.maskOpacity;
  }

  get mouseMode() {
    return this.state.mouseMode;
  }
}
