import { Mutations } from "vuex-smart-module";
import { groupListSchema, GroupState } from ".";
import { normalize } from "normalizr";
import { IGroup } from "./models";
import { IMember } from "@/modules/member/models";

export class GroupMutations extends Mutations<GroupState> {
  setActiveGroupId(id: number | null) {
    this.state.activeGroupId = id;
  }

  setEntities(payload: IGroup[]) {
    const normalizedData = normalize<IGroup>(payload, groupListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.groups ? Object.freeze(normalizedData.entities.groups) : {};
  }

  setEntity(payload: IGroup) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  addEntity(payload: IGroup) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  updateEntity(payload: IGroup) {
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  deleteEntity(id: number) {
    this.state.ids = this.state.ids.filter((item) => item !== id);
    const entities = { ...this.state.entities };
    delete entities[id];
    this.state.entities = Object.freeze(entities);
  }

  setMyMember(payload: IMember) {
    this.state.myMember = payload;
  }

  reset() {
    // acquire initial state
    const s = new GroupState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
