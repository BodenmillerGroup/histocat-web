import { Mutations } from "vuex-smart-module";
import { modelListSchema, ModelsState } from ".";
import { IModel } from "./models";
import { normalize } from "normalizr";

export class ModelsMutations extends Mutations<ModelsState> {
  setActiveModelId(id: number | null) {
    this.state.activeModelId = id;
  }

  setEntities(payload: IModel[]) {
    const normalizedData = normalize<IModel>(payload, modelListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.models ? normalizedData.entities.models : {};
  }

  setEntity(payload: IModel) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IModel) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IModel) {
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
    const s = new ModelsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
