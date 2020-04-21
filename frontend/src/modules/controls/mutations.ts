import { Mutations } from "vuex-smart-module";
import { ControlsState } from ".";

export class ControlsMutations extends Mutations<ControlsState> {
  setGraphInteractionMode(value: string) {
    this.state.graphInteractionMode = value;
  }
}
