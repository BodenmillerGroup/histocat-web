import { Getters } from "vuex-smart-module";
import { DatasetsState } from ".";

export class DatasetsGetters extends Getters<DatasetsState> {
  get datasets() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getDataset(id: number) {
    return this.state.entities[id];
  }

  get activeDatasetId() {
    return this.state.activeDatasetId;
  }

  get activeDataset() {
    return this.getters.activeDatasetId ? this.getters.getDataset(this.getters.activeDatasetId) : null;
  }

  get channels() {
    if (this.getters.activeDataset && this.getters.activeDataset.meta["channel_map"]) {
      const dict = this.getters.activeDataset.meta["channel_map"];
      return Object.keys(dict).sort((a, b) => {
        return dict[a] - dict[b];
      });
    }
    return [];
  }
}
