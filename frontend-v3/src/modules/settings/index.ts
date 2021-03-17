import create from "zustand";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";

type SettingsState = {
  channelsSettings: { [channelName: string]: IChannelSettings };
  filter: IImageFilter;
  legend: IImageLegend;
  scalebar: IImageScalebar;
  mask: IMaskSettings;
  mouseMode: "panZoom" | "lasso" | "rotate";

  setChannelsSettings(payload: any): void;
  setChannelColor(channelName: string, color: string): void;
  setChannelLevels(channelName: string, levels: { min: number; max: number }): void;
};

export const useSettingsStore = create<SettingsState>((set, get) => ({
  channelsSettings: {},
  filter: {
    apply: false,
    type: "gaussian",
    settings: {
      sigma: 1.0,
      kernel_size: 3,
    },
  },
  legend: {
    apply: false,
    fontScale: 12,
    showIntensity: false,
  },
  scalebar: {
    apply: false,
    settings: {
      scale: 1.0,
    },
  },
  mask: {
    mode: "raw",
  },
  mouseMode: "panZoom",

  setChannelsSettings(payload: any) {
    set({ channelsSettings: { ...get().channelsSettings, payload } });
  },

  setChannelColor(channelName: string, color: string) {
    const channelsSettings = get().channelsSettings;
    const existingLevels =
      channelsSettings[channelName] && channelsSettings[channelName].levels
        ? channelsSettings[channelName].levels
        : undefined;
    set({
      channelsSettings: {
        ...channelsSettings,
        [channelName]: { levels: existingLevels, color: color },
      },
    });
  },

  setChannelLevels(channelName: string, levels: { min: number; max: number }) {
    const channelsSettings = get().channelsSettings;
    const existingColor = channelsSettings[channelName] ? channelsSettings[channelName].color : "#ffffff";
    set({
      channelsSettings: {
        ...channelsSettings,
        [channelName]: { levels: levels, color: existingColor },
      },
    });
  },
}));
