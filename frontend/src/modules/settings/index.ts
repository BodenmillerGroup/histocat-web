import { Module } from "vuex-smart-module";
import { SettingsActions } from "./actions";
import { SettingsGetters } from "./getters";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar } from "./models";
import { SettingsMutations } from "./mutations";
import { DEFAULT_LAYOUT_UID } from "@/modules/ui/constants";

export class SettingsState {
  activeLayoutUid = DEFAULT_LAYOUT_UID; // Default layout uid

  channelsSettings: { [channelName: string]: IChannelSettings } = {};
  filter: IImageFilter = {
    apply: false,
    type: "gaussian",
    settings: {
      sigma: 1.0,
      kernel_size: 3,
    },
  };
  legend: IImageLegend = {
    apply: false,
    fontScale: 12,
    showIntensity: false,
  };
  scalebar: IImageScalebar = {
    apply: false,
    settings: {
      scale: 1.0,
    },
  };
}

export const settingsModule = new Module({
  namespaced: true,

  state: SettingsState,
  getters: SettingsGetters,
  mutations: SettingsMutations,
  actions: SettingsActions,
});
