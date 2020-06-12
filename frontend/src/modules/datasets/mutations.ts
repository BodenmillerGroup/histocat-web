import { WebSocketMessage } from "@/utils/WebSocketMessage";
import { Mutations } from "vuex-smart-module";
import { datasetListSchema, DatasetState } from ".";
import { IDataset } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_ACTIVE_DATASET_ID, SET_DATASETS } from "@/modules/datasets/events";
import { normalize } from "normalizr";

export class DatasetMutations extends Mutations<DatasetState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_DATASETS, (payload) => this.setEntities(payload));
    BroadcastManager.subscribe(SET_ACTIVE_DATASET_ID, (payload) => this.setActiveDatasetId(payload));
  }

  setActiveDatasetId(id: number | null) {
    this.state.activeDatasetId = id;
  }

  setEntities(payload: IDataset[]) {
    const normalizedData = normalize<IDataset>(payload, datasetListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.datasets ? normalizedData.entities.datasets : {};
  }

  setEntity(payload: IDataset) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IDataset) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IDataset) {
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  deleteEntity(id: number) {
    this.state.ids = this.state.ids.filter((item) => item !== id);
    const entities = { ...this.state.entities };
    delete entities[id];
    this.state.entities = entities;
  }

  updateDatasetTSNEOutput(message: WebSocketMessage) {
    const newState = { ...this.state.entities };
    const dataset = newState[message.payload.params.dataset_id];
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
      this.state.entities = newState;
    }
  }

  updateDatasetUMAPOutput(message: WebSocketMessage) {
    const newState = { ...this.state.entities };
    const dataset = newState[message.payload.params.dataset_id];
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
      this.state.entities = newState;
    }
  }

  updateDatasetPhenoGraphOutput(message: WebSocketMessage) {
    const newState = { ...this.state.entities };
    const dataset = newState[message.payload.params.dataset_id];
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
      this.state.entities = newState;
    }
  }

  reset() {
    // acquire initial state
    const s = new DatasetState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
