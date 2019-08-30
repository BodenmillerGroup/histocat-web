import { IImageSegmentationSettings, ImageResultType } from '@/modules/analysis/models';
import { Module } from 'vuex-smart-module';
import { SettingsActions } from './actions';
import { SettingsGetters } from './getters';
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from './models';
import { SettingsMutations } from './mutations';

export class SettingsState {
  channelsSettings: Map<number, IChannelSettings> = new Map<number, IChannelSettings>();
  metalColorMap: Map<string, string> = new Map<string, string>();
  filter: IImageFilter = {
    apply: false,
    type: 'gaussian',
    settings: {
      sigma: 1.0,
      kernel_size: 1,
    },
  };
  legend: IImageLegend = {
    apply: false,
    fontScale: 1.0,
    showIntensity: true,
  };
  scalebar: IImageScalebar = {
    apply: false,
    settings: {
      scale: 1.0,
    },
  };
  segmentationSettings: IImageSegmentationSettings = {
    algorithm: 'Otsu Hue',
    iterations: 1,
    kernel_size: 3,
    mask_color: '#00AAFF40',
    result_type: 'origin'
  };
  mask: IMaskSettings = {
    apply: false,
    location: undefined,
  }
}

export const settingsModule = new Module({
  namespaced: false,

  state: SettingsState,
  getters: SettingsGetters,
  mutations: SettingsMutations,
  actions: SettingsActions,
});
