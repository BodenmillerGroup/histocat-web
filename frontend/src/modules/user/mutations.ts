import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { UserState } from '.';
import { IUserProfile } from './models';

export const mutations = {
  setUsers(state: UserState, payload: IUserProfile[]) {
    state.users = payload;
  },
  setUser(state: UserState, payload: IUserProfile) {
    const users = state.users.filter((user: IUserProfile) => user.id !== payload.id);
    users.push(payload);
    state.users = users;
  },
};

const { commit } = getStoreAccessors<UserState, RootState>('');

export const commitSetUser = commit(mutations.setUser);
export const commitSetUsers = commit(mutations.setUsers);
