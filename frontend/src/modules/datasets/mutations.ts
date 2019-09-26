import { WebSocketMessage } from "@/utils/WebSocketMessage";
import { Mutations } from "vuex-smart-module";
import { DatasetState } from ".";
import { IDataset } from "./models";

export class DatasetMutations extends Mutations<DatasetState> {
  setDatasets(datasets: IDataset[]) {
    this.state.datasets = datasets;
  }

  setDataset(dataset: IDataset) {
    const items = this.state.datasets.filter(item => item.id !== dataset.id);
    items.push(dataset);
    this.state.datasets = items;
  }

  deleteDataset(id: number) {
    this.state.datasets = this.state.datasets.filter(item => item.id !== id);
  }

  setActiveDataset(dataset?: IDataset) {
    this.state.activeDataset = dataset;
  }

  updateDatasetTSNEOutput(message: WebSocketMessage) {
    const dataset = this.state.datasets.find(item => item.id === message.payload.params.dataset_id);
    if (dataset) {
      if (!dataset.output) {
        dataset.output = {
          tsne: {},
          umap: {},
          phenograph: {},
        };
      }
      if (!dataset.output.tsne) {
        dataset.output.tsne = {};
      }
      dataset.output.tsne[message.payload.name] = message.payload;
      this.state.activeDataset = Object.assign({}, dataset);
    }
  }

  updateDatasetUMAPOutput(message: WebSocketMessage) {
    const dataset = this.state.datasets.find(item => item.id === message.payload.params.dataset_id);
    if (dataset) {
      if (!dataset.output) {
        dataset.output = {
          tsne: {},
          umap: {},
          phenograph: {},
        };
      }
      if (!dataset.output.umap) {
        dataset.output.umap = {};
      }
      dataset.output.umap[message.payload.name] = message.payload;
      this.state.activeDataset = Object.assign({}, dataset);
    }
  }

  updateDatasetPhenoGraphOutput(message: WebSocketMessage) {
    const dataset = this.state.datasets.find(item => item.id === message.payload.params.dataset_id);
    if (dataset) {
      if (!dataset.output) {
        dataset.output = {
          tsne: {},
          umap: {},
          phenograph: {},
        };
      }
      if (!dataset.output.phenograph) {
        dataset.output.phenograph = {};
      }
      dataset.output.phenograph[message.payload.name] = message.payload;
      this.state.activeDataset = Object.assign({}, dataset);
    }
  }

  reset() {
    this.state.datasets = [];
    this.state.activeDataset = undefined;
  }
}
