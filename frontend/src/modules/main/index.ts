import { mutations } from './mutations';
import { getters } from './getters';
import { actions } from './actions';
import { MainState } from './state';
import { Module } from 'vuex';
import { RootState } from '@/store/state';

const defaultState: MainState = {
  isLoggedIn: null,
  token: '',
  logInError: false,
  userProfile: null,
  dashboardMiniDrawer: true,
  dashboardShowDrawer: true,
  notifications: [],
};

export const mainModule: Module<MainState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
