import { Module } from 'vuex-smart-module';
import { SettingsActions } from './actions';
import { SettingsGetters } from './getters';
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar } from './models';
import { SettingsMutations } from './mutations';

export class SettingsState {
  channelsSettings: Map<number, IChannelSettings> = new Map<number, IChannelSettings>();
  metalColorMap: Map<string, string> = new Map<string, string>();
  filter: IImageFilter = {
    apply: false,
    type: 'gaussian',
    settings: {
      sigma: 1.0
    },
  };
  legend: IImageLegend = {
    apply: false,
    fontScale: 1.0,
  };
  scalebar: IImageScalebar = {
    apply: false,
    settings: {
      scale: 1.0
    },
  };
}

export const settingsModule = new Module({
  namespaced: false,

  state: SettingsState,
  getters: SettingsGetters,
  mutations: SettingsMutations,
  actions: SettingsActions,
});
