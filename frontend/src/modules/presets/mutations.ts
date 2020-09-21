import { Mutations } from "vuex-smart-module";
import { presetListSchema, PresetsState } from ".";
import { IPreset } from "./models";
import { normalize } from "normalizr";

export class PresetsMutations extends Mutations<PresetsState> {
  setActivePresetId(id: number | null) {
    this.state.activePresetId = id;
  }

  setEntities(payload: IPreset[]) {
    const normalizedData = normalize<IPreset>(payload, presetListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.presets ? normalizedData.entities.presets : {};
  }

  setEntity(payload: IPreset) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IPreset) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IPreset) {
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
    const s = new PresetsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
