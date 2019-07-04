import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { SettingsState } from '.';

export const getters = {
  channelSettings: (state: SettingsState) => (id: number) => {
    return state.channelsSettings.get(id);
  },
  metalColorMap: (state: SettingsState) => {
    return state.metalColorMap;
  },
};

const { read } = getStoreAccessors<SettingsState, RootState>('');

export const readChannelSettings = read(getters.channelSettings);
export const readMetalColorMap = read(getters.metalColorMap);