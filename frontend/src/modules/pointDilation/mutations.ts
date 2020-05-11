import { Mutations } from "vuex-smart-module";
import { PointDilationState } from ".";

export class PointDilationMutations extends Mutations<PointDilationState> {
  setMetadataField(value: string) {
    this.state.metadataField = value;
  }

  setCategoryField(value: string) {
    this.state.categoryField = value;
  }
}
