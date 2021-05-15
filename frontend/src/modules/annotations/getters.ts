import { Getters } from "vuex-smart-module";
import { AnnotationsState } from ".";

export class AnnotationsGetters extends Getters<AnnotationsState> {
  get cellClasses() {
    return this.state.cellClasses;
  }

  get annotations() {
    return this.state.annotations;
  }
}
