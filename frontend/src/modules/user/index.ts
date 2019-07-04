import { RootState } from '@/store';
import { Module } from 'vuex';
import { actions } from './actions';
import { getters } from './getters';
import { IUserProfile } from './models';
import { mutations } from './mutations';

export interface UserState {
  users: IUserProfile[];
}

const defaultState: UserState = {
  users: [],
};

export const userModule: Module<UserState, RootState> = {
  state: defaultState,
  mutations,
  actions,
  getters,
};
