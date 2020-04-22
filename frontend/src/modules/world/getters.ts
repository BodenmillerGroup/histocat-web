import { Getters } from "vuex-smart-module";
import { WorldState } from ".";

export class WorldGetters extends Getters<WorldState> {
  get world() {
    return this.state.world;
  }
}
