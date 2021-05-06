import { Mutations } from "vuex-smart-module";
import { PROJECT_LAYOUTS_STORAGE_KEY, UiState } from ".";
import { ILayout, IResponsive } from "./models";
import { GoldenLayout, LayoutConfig } from "golden-layout";
import { DEFAULT_LAYOUTS } from "@/modules/ui/defaultLayouts";

export class UiMutations extends Mutations<UiState> {
  setGoldenLayout(value: GoldenLayout | null) {
    this.state.goldenLayout = value;
  }

  addLayout(name: string) {
    if (this.state.goldenLayout) {
      const layout: ILayout = {
        name: name,
        config: LayoutConfig.fromResolved(this.state.goldenLayout.saveLayout()),
        isDefault: false,
      };
      const layouts = this.state.layouts.concat(layout);
      this.state.layouts = layouts;
      localStorage.setItem(PROJECT_LAYOUTS_STORAGE_KEY, JSON.stringify(layouts));
    }
  }

  resetLayouts() {
    const layouts = DEFAULT_LAYOUTS;
    this.state.layouts = layouts;
    this.state.activeLayout = layouts[0];
    localStorage.setItem(PROJECT_LAYOUTS_STORAGE_KEY, JSON.stringify(layouts));
    if (this.state.goldenLayout) {
      this.state.goldenLayout.loadLayout(this.state.activeLayout.config);
    }
  }

  loadLayout(name: string) {
    if (this.state.goldenLayout) {
      const layout = this.state.layouts.find((item) => item.name === name);
      if (layout) {
        this.state.activeLayout = layout;
        this.state.goldenLayout.loadLayout(layout.config);
      }
    }
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
