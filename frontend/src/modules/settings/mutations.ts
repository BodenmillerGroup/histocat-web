import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IImageFilter, IImageLegend, IImageScalebar } from "./models";

export class SettingsMutations extends Mutations<SettingsState> {
  setActiveLayoutUid(uid: string) {
    this.state.activeLayoutUid = uid;
  }

  setChannelsSettings(payload) {
    this.state.channelsSettings = { ...this.state.channelsSettings, ...payload };
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

  reset() {
    // acquire initial state
    const s = new SettingsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
