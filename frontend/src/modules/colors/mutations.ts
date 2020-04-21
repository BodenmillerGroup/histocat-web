import { Mutations } from "vuex-smart-module";
import {ColorsState} from ".";

export class ColorsMutations extends Mutations<ColorsState> {
  setColorMode(value) {
    this.state.colorMode = value;
  }

  setColorAccessor(value) {
    this.state.colorAccessor = value;
  }

  setRgb(value) {
    this.state.rgb = value;
  }

  setScale(value) {
    this.state.scale = value;
  }
}
