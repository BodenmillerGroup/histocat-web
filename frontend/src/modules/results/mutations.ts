import { Mutations } from "vuex-smart-module";
import { resultListSchema, ResultsState } from ".";
import { IResult } from "./models";
import { normalize } from "normalizr";

export class ResultsMutations extends Mutations<ResultsState> {
  setActiveResultId(id: number | null) {
    this.state.activeResultId = id;
  }

  setEntities(payload: IResult[]) {
    const normalizedData = normalize<IResult>(payload, resultListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.results ? normalizedData.entities.results : {};
  }

  setEntity(payload: IResult) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IResult) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IResult) {
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
    const s = new ResultsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
