import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_CHANNEL_SETTINGS,
  SET_FILTER,
  SET_LEGEND,
  SET_MASK_SETTINGS,
  SET_SCALEBAR,
} from "@/modules/settings/events";

export class SettingsMutations extends Mutations<SettingsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_CHANNEL_SETTINGS, (payload) => this.setChannelSettings(payload));
    BroadcastManager.subscribe(SET_FILTER, (payload) => this.setFilter(payload));
    BroadcastManager.subscribe(SET_LEGEND, (payload) => this.setLegend(payload));
    BroadcastManager.subscribe(SET_SCALEBAR, (payload) => this.setScalebar(payload));
    BroadcastManager.subscribe(SET_MASK_SETTINGS, (payload) => this.setMaskSettings(payload));
  }

  setChannelSettings(payload: { channelName: string; settings: IChannelSettings }) {
    this.state.channelsSettings = { ...this.state.channelsSettings, [payload.channelName]: payload.settings };
  }

  setChannelColor(payload: { channelName: string; color: string }) {
    const existingLevels =
      this.state.channelsSettings[payload.channelName] && this.state.channelsSettings[payload.channelName].levels
        ? this.state.channelsSettings[payload.channelName].levels
        : undefined;
    this.state.channelsSettings = {
      ...this.state.channelsSettings,
      [payload.channelName]: { levels: existingLevels, color: payload.color },
    };
  }

  setChannelLevels(payload: { channelName: string; levels: { min: number; max: number } }) {
    const existingColor = this.state.channelsSettings[payload.channelName]
      ? this.state.channelsSettings[payload.channelName].color
      : "#ffffff";
    this.state.channelsSettings = {
      ...this.state.channelsSettings,
      [payload.channelName]: { levels: payload.levels, color: existingColor },
    };
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

  setMouseMode(mode: "panZoom" | "lasso" | "rotate") {
    this.state.mouseMode = mode;
  }

  reset() {
    // acquire initial state
    const s = new SettingsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
