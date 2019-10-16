import { IImageSegmentationSettings } from "@/modules/analysis/models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";

export class SettingsMutations extends Mutations<SettingsState> {
  resetSettings(payload: { suppressBroadcast?: boolean }) {
    this.state.channelsSettings = new Map<number, IChannelSettings>();
    this.state.metalColorMap = new Map<string, string>();
    this.state.filter = {
      apply: false,
      type: "gaussian",
      settings: {
        sigma: 1.0,
        kernel_size: 3
      }
    };
    this.state.legend = {
      apply: false,
      fontScale: 1.0,
      showIntensity: false
    };
    this.state.scalebar = {
      apply: false,
      settings: {
        scale: 1.0
      }
    };
    this.state.segmentationSettings = {
      algorithm: "Otsu Hue",
      iterations: 1,
      kernel_size: 3,
      mask_color: "#00AAFF40",
      result_type: "origin"
    };
    this.state.mask = {
      apply: false,
      location: undefined
    };
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.resetSettings.name,
        payload: payload
      });
    }
  }

  setChannelSettings(payload: IChannelSettings) {
    this.state.channelsSettings.set(payload.id, payload);
    this.state.channelsSettings = new Map(this.state.channelsSettings);
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setChannelSettings.name,
        payload: payload
      });
    }
  }

  setMetalColor(payload: { metal: string; color: string; suppressBroadcast?: boolean }) {
    this.state.metalColorMap.set(payload.metal, payload.color);
    this.state.metalColorMap = new Map(this.state.metalColorMap);
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setMetalColor.name,
        payload: payload
      });
    }
  }

  setFilter(payload: IImageFilter) {
    this.state.filter = payload;
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setFilter.name,
        payload: payload
      });
    }
  }

  setLegend(payload: IImageLegend) {
    this.state.legend = payload;
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setLegend.name,
        payload: payload
      });
    }
  }

  setScalebar(payload: IImageScalebar) {
    this.state.scalebar = payload;
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setScalebar.name,
        payload: payload
      });
    }
  }

  setSegmentationSettings(payload: IImageSegmentationSettings) {
    this.state.segmentationSettings = payload;
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setSegmentationSettings.name,
        payload: payload
      });
    }
  }

  setMaskSettings(payload: IMaskSettings) {
    this.state.mask = payload;
    if (!payload.suppressBroadcast) {
      BroadcastManager.postMessage({
        method: this.setMaskSettings.name,
        payload: payload
      });
    }
  }
}
