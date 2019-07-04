import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { SettingsState } from '.';
import { IChannelSettings } from './models';

export const mutations = {
  setChannelSettings(state: SettingsState, channelSettings: IChannelSettings) {
    state.channelsSettings.set(channelSettings.id, channelSettings);
  },
  setMetalColor(state: SettingsState, payload: { metal: string, color: string }) {
    state.metalColorMap.set(payload.metal, payload.color);
    state.metalColorMap = new Map(state.metalColorMap);
  },
};

const { commit } = getStoreAccessors<SettingsState, RootState>('');

export const commitSetChannelSettings = commit(mutations.setChannelSettings);
export const commitSetMetalColor = commit(mutations.setMetalColor);
