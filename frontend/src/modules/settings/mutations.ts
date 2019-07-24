import { Mutations } from 'vuex-smart-module';
import { SettingsState } from '.';
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar } from './models';


export class SettingsMutations extends Mutations<SettingsState> {
  resetSettings() {
    this.state.channelsSettings = new Map<number, IChannelSettings>();
    this.state.metalColorMap = new Map<string, string>();
    this.state.filter = {
      apply: false,
      type: 'gaussian',
      settings: {},
    };
    this.state.legend = {
      apply: false,
      settings: {},
    };
  }

  setChannelSettings(channelSettings: IChannelSettings) {
    this.state.channelsSettings.set(channelSettings.id, channelSettings);
  }

  setMetalColor(payload: { metal: string, color: string }) {
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
}
