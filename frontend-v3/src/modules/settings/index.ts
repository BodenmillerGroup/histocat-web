import create from "zustand";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";

type SettingsState = {
  channelsSettings: { [channelName: string]: IChannelSettings };
  filter: IImageFilter;
  legend: IImageLegend;
  scalebar: IImageScalebar;
  mask: IMaskSettings;
  mouseMode: "panZoom" | "lasso" | "rotate";
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
}));
