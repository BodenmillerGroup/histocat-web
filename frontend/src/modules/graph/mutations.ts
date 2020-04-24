import { Mutations } from "vuex-smart-module";
import { GraphState } from ".";
import Dataframe from "@/cellxgene/util/dataframe/dataframe";

export class GraphMutations extends Mutations<GraphState> {
  setSchema(value: any) {
    this.state.schema = value;
  }

  setObsAnnotations(value: Dataframe) {
    this.state.obsAnnotations = value;
  }

  setVarAnnotations(value: Dataframe) {
    this.state.varAnnotations = value;
  }
}
