import { Getters } from "vuex-smart-module";
import { AnnotationsState } from ".";

export class AnnotationsGetters extends Getters<AnnotationsState> {
  get classes() {
    return this.state.classes;
  }

  get annotations() {
    return this.state.annotations;
  }
}
