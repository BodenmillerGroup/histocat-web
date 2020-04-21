import { Mutations } from "vuex-smart-module";
import { LayoutChoiceState } from ".";

export class LayoutChoiceMutations extends Mutations<LayoutChoiceState> {
  setAvailable(value: string[]) {
    this.state.available = value;
  }

  setCurrent(value?: string) {
    this.state.current = value;
  }

  setCurrentDimNames(value: string[]) {
    this.state.currentDimNames = value;
  }
}
