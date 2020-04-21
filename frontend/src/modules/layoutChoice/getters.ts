import { Getters } from "vuex-smart-module";
import { LayoutChoiceState } from ".";

export class LayoutChoiceGetters extends Getters<LayoutChoiceState> {
  get available() {
    return this.state.available;
  }

  get current() {
    return this.state.current;
  }

  get currentDimNames() {
    return this.state.currentDimNames;
  }
}
