import { Getters } from 'vuex-smart-module';
import { UserState } from '.';

export class UserGetters extends Getters<UserState> {
  get adminUsers() {
    return this.state.users;
  }

  adminOneUser(userId: number) {
    const filteredUsers = this.state.users.filter((user) => user.id === userId);
    if (filteredUsers.length > 0) {
      return { ...filteredUsers[0] };
    }
  }
}
