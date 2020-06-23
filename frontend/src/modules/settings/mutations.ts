import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_PRESET,
  SET_CHANNEL_SETTINGS,
  SET_FILTER,
  SET_LEGEND,
  SET_MASK_SETTINGS,
  SET_METAL_COLOR,
  SET_SCALEBAR,
  SET_SHARED_CHANNEL_SETTINGS,
} from "@/modules/settings/events";
import { IPreset } from "@/modules/presets/models";

export class SettingsMutations extends Mutations<SettingsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_SHARED_CHANNEL_SETTINGS, (payload) => this.setSharedChannelSettings(payload));
    BroadcastManager.subscribe(SET_CHANNEL_SETTINGS, (payload) => this.setChannelSettings(payload));
    BroadcastManager.subscribe(SET_METAL_COLOR, (payload) => this.setMetalColor(payload));
    BroadcastManager.subscribe(SET_FILTER, (payload) => this.setFilter(payload));
    BroadcastManager.subscribe(SET_LEGEND, (payload) => this.setLegend(payload));
    BroadcastManager.subscribe(SET_SCALEBAR, (payload) => this.setScalebar(payload));
    BroadcastManager.subscribe(SET_MASK_SETTINGS, (payload) => this.setMaskSettings(payload));
    BroadcastManager.subscribe(SET_PRESET, (payload) => this.setPreset(payload));
  }

  setSharedChannelSettings(payload: IChannelSettings[]) {
    let newState = { ...this.state.channelSettings };
    for (const item of payload) {
      let acquisitionChannelSettings = this.state.channelSettings[item.acquisitionId];
      if (!acquisitionChannelSettings) {
        acquisitionChannelSettings = {};
      }
      acquisitionChannelSettings[item.name] = item;
      newState[item.acquisitionId] = acquisitionChannelSettings;
    }
    this.state.channelSettings = newState;
  }

  setChannelSettings(payload: IChannelSettings) {
    let acquisitionChannelSettings = this.state.channelSettings[payload.acquisitionId];
    if (!acquisitionChannelSettings) {
      acquisitionChannelSettings = {};
    }
    acquisitionChannelSettings[payload.name] = payload;
    this.state.channelSettings = { ...this.state.channelSettings, [payload.acquisitionId]: acquisitionChannelSettings };
  }

  setMetalColor(payload: { metal: string; color: string }) {
    this.state.colorMap = { ...this.state.colorMap, [payload.metal]: payload.color };
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

  setMaskSettings(payload: IMaskSettings) {
    this.state.mask = payload;
  }

  setPreset(preset: IPreset) {
    this.state.channelSettings = preset.data["channelSettings"];
    this.state.colorMap = preset.data["colorMap"];
  }

  reset() {
    // acquire initial state
    const s = new SettingsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
