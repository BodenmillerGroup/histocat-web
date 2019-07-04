import { IUserProfile } from '@/modules/user/models';
import { RootState } from '@/store';
import { Module } from 'vuex';
import { actions } from './actions';
import { getters } from './getters';
import { AppNotification } from './models';
import { mutations } from './mutations';

export interface MainState {
  token: string;
  isLoggedIn: boolean | null;
  logInError: boolean;
  userProfile: IUserProfile | null;
  dashboardMiniDrawer: boolean;
  dashboardShowDrawer: boolean;
  notifications: AppNotification[];
}

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
