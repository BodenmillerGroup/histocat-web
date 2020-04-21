import { Getters } from "vuex-smart-module";
import { PointDilationState } from ".";

export class PointDilationGetters extends Getters<PointDilationState> {
  get metadataField() {
    return this.state.metadataField;
  }

  get categoryField() {
    return this.state.categoryField;
  }
}
