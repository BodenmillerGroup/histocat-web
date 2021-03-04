import { Mutations } from "vuex-smart-module";
import { memberListSchema, MemberState } from ".";
import { normalize } from "normalizr";
import { IMember } from "./models";

export class MemberMutations extends Mutations<MemberState> {
  setEntities(payload: IMember[]) {
    const normalizedData = normalize<IMember>(payload, memberListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.members ? normalizedData.entities.members : {};
  }

  setEntity(payload: IMember) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IMember) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IMember) {
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
    const s = new MemberState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
