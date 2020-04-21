import { Getters } from "vuex-smart-module";
import { ColorsState } from ".";

export class ColorsGetters extends Getters<ColorsState> {
  get colorMode() {
    return this.state.colorMode;
  }

  get colorAccessor() {
    return this.state.colorAccessor;
  }

  get rgb() {
    return this.state.rgb;
  }

  get scale() {
    return this.state.scale;
  }
}
