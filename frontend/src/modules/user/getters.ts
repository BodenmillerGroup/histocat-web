import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { UserState } from '.';

export const getters = {
  adminUsers: (state: UserState) => state.users,
  adminOneUser: (state: UserState) => (userId: number) => {
    const filteredUsers = state.users.filter((user) => user.id === userId);
    if (filteredUsers.length > 0) {
      return { ...filteredUsers[0] };
    }
  },
};

const { read } = getStoreAccessors<UserState, RootState>('');

export const readAdminOneUser = read(getters.adminOneUser);
export const readAdminUsers = read(getters.adminUsers);
