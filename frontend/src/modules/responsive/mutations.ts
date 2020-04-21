import { Mutations } from "vuex-smart-module";
import { IResponsive } from "./models";
import { ResponsiveState } from ".";

export class ResponsiveMutations extends Mutations<ResponsiveState> {
  setResponsive(value: IResponsive) {
    this.state.responsive = value;
  }
}
