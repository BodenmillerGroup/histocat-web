import { Getters } from "vuex-smart-module";
import { UniverseState } from ".";

export class UniverseGetters extends Getters<UniverseState> {
  get universe() {
    return this.state.universe;
  }
}
