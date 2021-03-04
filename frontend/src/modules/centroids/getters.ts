import { Getters } from "vuex-smart-module";
import { CentroidsState } from ".";

export class CentroidsGetters extends Getters<CentroidsState> {
  get centroids() {
    return this.state.centroids;
  }
}
