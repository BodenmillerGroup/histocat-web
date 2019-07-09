import { Module } from 'vuex-smart-module';
import { SettingsActions } from './actions';
import { SettingsGetters } from './getters';
import { IChannelSettings } from './models';
import { SettingsMutations } from './mutations';

export class SettingsState {
  channelsSettings: Map<number, IChannelSettings> = new Map<number, IChannelSettings>();
  metalColorMap: Map<string, string> = new Map<string, string>();
}

export const settingsModule = new Module({
  namespaced: false,

  state: SettingsState,
  getters: SettingsGetters,
  mutations: SettingsMutations,
  actions: SettingsActions,
});