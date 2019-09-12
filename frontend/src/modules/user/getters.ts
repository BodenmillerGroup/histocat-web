import { Getters } from "vuex-smart-module";
import { UserState } from ".";

export class UserGetters extends Getters<UserState> {
  get users() {
    return this.state.users;
  }

  getUser(userId: number) {
    const filteredUsers = this.users.filter(user => user.id === userId);
    if (filteredUsers.length > 0) {
      return { ...filteredUsers[0] };
    }
  }
}
