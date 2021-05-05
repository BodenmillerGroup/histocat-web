import { Mutations } from "vuex-smart-module";
import { UiState } from ".";
import { ILayout, IResponsive } from "./models";

export class UiMutations extends Mutations<UiState> {
  setActiveLayout(value: ILayout) {
    this.state.activeLayout = value;
  }

  setResponsive(value: IResponsive) {
    this.state.responsive = value;
  }

  setProcessing(payload: boolean) {
    this.state.processing = payload;
  }

  setProcessingProgress(payload: number) {
    this.state.processingProgress = payload;
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
