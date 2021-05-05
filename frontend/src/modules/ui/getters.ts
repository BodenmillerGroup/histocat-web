import { Getters } from "vuex-smart-module";
import { UiState } from ".";

export class UiGetters extends Getters<UiState> {
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
