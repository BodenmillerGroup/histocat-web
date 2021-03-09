import { Mutations } from "vuex-smart-module";
import { datasetListSchema, DatasetsState } from ".";
import { IDataset } from "./models";
import { normalize } from "normalizr";

export class DatasetsMutations extends Mutations<DatasetsState> {
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

  reset() {
    // acquire initial state
    const s = new DatasetsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
