import { IChannelSettings } from '@/modules/settings/models';
import { RootState } from '@/store';
import { Module } from 'vuex';
import { actions } from './actions';
import { getters } from './getters';
import { mutations } from './mutations';

export interface SettingsState {
  channelsSettings: Map<number, IChannelSettings>;
  metalColorMap: Map<string, string>;
}

const defaultState: SettingsState = {
  channelsSettings: new Map<number, IChannelSettings>(),
  metalColorMap: new Map<string, string>(),
};

export const settingsModule: Module<SettingsState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
