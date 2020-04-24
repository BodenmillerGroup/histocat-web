import { Getters } from "vuex-smart-module";
import { GraphState } from ".";

export class GraphGetters extends Getters<GraphState> {
  get schema() {
    return this.state.schema;
  }

  get obsAnnotations() {
    return this.state.obsAnnotations;
  }

  get varAnnotations() {
    return this.state.varAnnotations;
  }
}
