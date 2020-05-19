import { IImageSegmentationSettings } from "@/modules/analysis/models";
import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  RESET_SETTINGS,
  SET_CHANNEL_SETTINGS,
  SET_FILTER,
  SET_LEGEND,
  SET_MASK_SETTINGS,
  SET_METAL_COLOR,
  SET_SCALEBAR,
  SET_SEGMENTATION_SETTINGS, SET_SHARED_CHANNEL_SETTINGS,
} from "@/modules/settings/events";

export class SettingsMutations extends Mutations<SettingsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_SHARED_CHANNEL_SETTINGS, (payload) => this.setSharedChannelSettings(payload));
    BroadcastManager.subscribe(SET_CHANNEL_SETTINGS, (payload) => this.setChannelSettings(payload));
    BroadcastManager.subscribe(SET_METAL_COLOR, (payload) => this.setMetalColor(payload));
    BroadcastManager.subscribe(SET_FILTER, (payload) => this.setFilter(payload));
    BroadcastManager.subscribe(SET_LEGEND, (payload) => this.setLegend(payload));
    BroadcastManager.subscribe(SET_SCALEBAR, (payload) => this.setScalebar(payload));
    BroadcastManager.subscribe(SET_SEGMENTATION_SETTINGS, (payload) => this.setSegmentationSettings(payload));
    BroadcastManager.subscribe(SET_MASK_SETTINGS, (payload) => this.setMaskSettings(payload));
    BroadcastManager.subscribe(RESET_SETTINGS, (payload) => this.resetSettings());
  }

  setSharedChannelSettings(payload: IChannelSettings[]) {
    for (const item of payload) {
      let acquisitionChannelSettings = this.state.channelSettings.get(item.acquisitionId);
      if (!acquisitionChannelSettings) {
        acquisitionChannelSettings = new Map<string, IChannelSettings>();
      }
      acquisitionChannelSettings.set(item.name, item);
      this.state.channelSettings.set(item.acquisitionId, acquisitionChannelSettings);
    }
    this.state.channelSettings = new Map(this.state.channelSettings);
  }

  setChannelSettings(payload: IChannelSettings) {
    let acquisitionChannelSettings = this.state.channelSettings.get(payload.acquisitionId);
    if (!acquisitionChannelSettings) {
      acquisitionChannelSettings = new Map<string, IChannelSettings>();
    }
    acquisitionChannelSettings.set(payload.name, payload);
    this.state.channelSettings.set(payload.acquisitionId, acquisitionChannelSettings);
    this.state.channelSettings = new Map(this.state.channelSettings);
  }

  setMetalColor(payload: { metal: string; color: string; suppressBroadcast?: boolean }) {
    this.state.metalColorMap.set(payload.metal, payload.color);
    this.state.metalColorMap = new Map(this.state.metalColorMap);
  }

  setFilter(payload: IImageFilter) {
    this.state.filter = payload;
  }

  setLegend(payload: IImageLegend) {
    this.state.legend = payload;
  }

  setScalebar(payload: IImageScalebar) {
    this.state.scalebar = payload;
  }

  setSegmentationSettings(payload: IImageSegmentationSettings) {
    this.state.segmentationSettings = payload;
  }

  setMaskSettings(payload: IMaskSettings) {
    this.state.mask = payload;
  }

  resetSettings() {
    this.state.channelSettings = new Map<number, Map<string, IChannelSettings>>();
    this.state.metalColorMap = new Map<string, string>();
    this.state.filter = {
      apply: false,
      type: "gaussian",
      settings: {
        sigma: 1.0,
        kernel_size: 3,
      },
    };
    this.state.legend = {
      apply: false,
      fontScale: 1.0,
      showIntensity: false,
    };
    this.state.scalebar = {
      apply: false,
      settings: {
        scale: 1.0,
      },
    };
    this.state.segmentationSettings = {
      algorithm: "Otsu Hue",
      iterations: 1,
      kernel_size: 3,
      mask_color: "#00AAFF40",
      result_type: "origin",
    };
    this.state.mask = {
      apply: false,
      location: undefined,
    };
  }
}
