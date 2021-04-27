import { Getters } from "vuex-smart-module";
import { DatasetsState } from ".";

export class DatasetsGetters extends Getters<DatasetsState> {
  get datasets() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  get activeDatasetId() {
    return this.state.activeDatasetId;
  }

  get activeDataset() {
    return this.state.activeDatasetId ? this.state.entities[this.state.activeDatasetId] : null;
  }

  get channels() {
    return this.getters.activeDataset ? this.getters.activeDataset.channels : [];
  }
}
