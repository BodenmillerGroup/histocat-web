import { IImageSegmentationSettings } from "@/modules/analysis/models";
import { Mutations } from "vuex-smart-module";
import { SettingsState } from ".";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";

export class SettingsMutations extends Mutations<SettingsState> {
  resetSettings() {
    this.state.channelsSettings = new Map<number, IChannelSettings>();
    this.state.metalColorMap = new Map<string, string>();
    this.state.filter = {
      apply: false,
      type: "gaussian",
      settings: {
        sigma: 1.0,
        kernel_size: 1
      }
    };
    this.state.legend = {
      apply: false,
      fontScale: 1.0,
      showIntensity: true
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
  }

  setChannelSettings(channelSettings: IChannelSettings) {
    this.state.channelsSettings.set(channelSettings.id, channelSettings);
    this.state.channelsSettings = new Map(this.state.channelsSettings);
  }

  setMetalColor(payload: { metal: string; color: string }) {
    this.state.metalColorMap.set(payload.metal, payload.color);
    this.state.metalColorMap = new Map(this.state.metalColorMap);
  }

  setFilter(filter: IImageFilter) {
    this.state.filter = filter;
  }

  setLegend(legend: IImageLegend) {
    this.state.legend = legend;
  }

  setScalebar(scalebar: IImageScalebar) {
    this.state.scalebar = scalebar;
  }

  setSegmentationSettings(settings: IImageSegmentationSettings) {
    this.state.segmentationSettings = settings;
  }

  setMaskSettings(settings: IMaskSettings) {
    this.state.mask = settings;
  }
}
