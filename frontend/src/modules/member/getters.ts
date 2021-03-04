import { Getters } from "vuex-smart-module";
import { MemberState } from ".";

export class MemberGetters extends Getters<MemberState> {
  get members() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getMember(id: number) {
    return this.state.entities[id];
  }
}
