import { Mutations } from "vuex-smart-module";
import { IResponsive } from "./models";
import { ResponsiveState } from ".";

export class ResponsiveMutations extends Mutations<ResponsiveState> {
  setResponsive(value: IResponsive) {
    this.state.responsive = value;
  }

  reset() {
    // acquire initial state
    const s = new ResponsiveState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
