import { Mutations } from "vuex-smart-module";
import { gateListSchema, GateState } from ".";
import { IGate } from "./models";
import { normalize } from "normalizr";

export class GateMutations extends Mutations<GateState> {
  setActiveGateId(id: number | null) {
    this.state.activeGateId = id;
  }

  setEntities(payload: IGate[]) {
    const normalizedData = normalize<IGate>(payload, gateListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.gates ? normalizedData.entities.gates : {};
  }

  setEntity(payload: IGate) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IGate) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IGate) {
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
    const s = new GateState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
