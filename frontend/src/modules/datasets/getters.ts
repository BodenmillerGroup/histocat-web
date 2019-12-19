import { Getters } from "vuex-smart-module";
import { DatasetState } from ".";

export class DatasetGetters extends Getters<DatasetState> {
  get datasets() {
    return this.state.datasets;
  }

  get activeDataset() {
    return this.state.activeDataset;
  }

  get channels() {
    return this.getters.activeDataset && this.getters.activeDataset.input["channel_map"]
      ? Object.keys(this.getters.activeDataset.input["channel_map"])
      : [];
  }

  get heatmaps() {
    if (!this.getters.activeDataset || !this.getters.activeDataset.input["neighbors_columns"]) {
      return [];
    }
    const channelItems = this.getters.channels.map(item => {
      return {
        type: "channel",
        label: item
      };
    });
    const neighborItems = this.getters.activeDataset.input["neighbors_columns"].map(item => {
      return {
        type: "neighbor",
        label: item.substring(10, item.length)
      };
    });
    return channelItems.concat(neighborItems);
  }

  get tsneResults() {
    return this.getters.activeDataset && this.getters.activeDataset.output && this.getters.activeDataset.output["tsne"]
      ? Object.values(this.getters.activeDataset.output["tsne"])
      : [];
  }

  get umapResults() {
    return this.getters.activeDataset && this.getters.activeDataset.output && this.getters.activeDataset.output["umap"]
      ? Object.values(this.getters.activeDataset.output["umap"])
      : [];
  }

  get phenographResults() {
    return this.getters.activeDataset && this.getters.activeDataset.output && this.getters.activeDataset.output["phenograph"]
      ? Object.values(this.getters.activeDataset.output["phenograph"])
      : [];
  }
}
