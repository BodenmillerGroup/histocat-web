import { Getters } from "vuex-smart-module";
import { ResponsiveState } from ".";

export class ResponsiveGetters extends Getters<ResponsiveState> {
  get responsive() {
    return this.state.responsive;
  }
}
