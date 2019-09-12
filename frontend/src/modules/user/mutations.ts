import { Mutations } from "vuex-smart-module";
import { UserState } from ".";
import { IUserProfile } from "./models";

export class UserMutations extends Mutations<UserState> {
  setUsers(payload: IUserProfile[]) {
    this.state.users = payload;
  }

  setUser(payload: IUserProfile) {
    const users = this.state.users.filter((user: IUserProfile) => user.id !== payload.id);
    users.push(payload);
    this.state.users = users;
  }
}
