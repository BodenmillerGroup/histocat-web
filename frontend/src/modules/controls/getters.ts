import { Getters } from "vuex-smart-module";
import { ControlsState } from ".";

export class ControlsGetters extends Getters<ControlsState> {
  get graphInteractionMode() {
    return this.state.graphInteractionMode;
  }
}
