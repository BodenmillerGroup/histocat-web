import { IImageSegmentationSettings } from "@/modules/analysis/models";
import { Module } from "vuex-smart-module";
import { SettingsActions } from "./actions";
import { SettingsGetters } from "./getters";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";
import { SettingsMutations } from "./mutations";

export class SettingsState {
  channelSettings: { [acquisitionId: number]: { [channelName: string]: IChannelSettings } } = {};
  colorMap: { [markerName: string]: string } = {};
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
  segmentationSettings: IImageSegmentationSettings = {
    algorithm: "Otsu Hue",
    iterations: 1,
    kernel_size: 3,
    mask_color: "#00AAFF40",
    result_type: "origin",
  };
  mask: IMaskSettings = {
    apply: false,
    location: undefined,
  };
}

export const settingsModule = new Module({
  namespaced: true,

  state: SettingsState,
  getters: SettingsGetters,
  mutations: SettingsMutations,
  actions: SettingsActions,
});
