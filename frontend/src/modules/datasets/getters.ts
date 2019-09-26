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
    return this.activeDataset && this.activeDataset.input["channel_map"]
      ? Object.keys(this.activeDataset.input["channel_map"])
      : [];
  }

  get heatmaps() {
    if (!this.activeDataset || !this.activeDataset.input["neighbors_columns"]) {
      return [];
    }
    const channelItems = this.channels.map(item => {
      return {
        type: "channel",
        label: item
      };
    });
    const neighborItems = this.activeDataset.input["neighbors_columns"].map(item => {
      return {
        type: "neighbor",
        label: item.substring(10, item.length)
      };
    });
    return channelItems.concat(neighborItems);
  }

  get tsneResults() {
    return this.activeDataset && this.activeDataset.output && this.activeDataset.output["tsne"]
      ? Object.values(this.activeDataset.output["tsne"])
      : [];
  }

  get umapResults() {
    return this.activeDataset && this.activeDataset.output && this.activeDataset.output["umap"]
      ? Object.values(this.activeDataset.output["umap"])
      : [];
  }

  get phenographResults() {
    return this.activeDataset && this.activeDataset.output && this.activeDataset.output["phenograph"]
      ? Object.values(this.activeDataset.output["phenograph"])
      : [];
  }
}
