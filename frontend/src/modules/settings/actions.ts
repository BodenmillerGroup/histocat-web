import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { ActionContext } from 'vuex';
import { SettingsState } from '.';

type SettingsContext = ActionContext<SettingsState, RootState>;

export const actions = {};

const { dispatch } = getStoreAccessors<SettingsState, RootState>('');
