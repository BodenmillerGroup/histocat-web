import { mutations } from './mutations';
import { getters } from './getters';
import { actions } from './actions';
import { UserState } from './state';
import { Module } from 'vuex';
import { RootState } from '@/store/state';

const defaultState: UserState = {
  users: [],
};

export const userModule: Module<UserState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
