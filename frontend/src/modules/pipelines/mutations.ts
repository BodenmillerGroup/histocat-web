import { Mutations } from "vuex-smart-module";
import { pipelineListSchema, PipelinesState } from ".";
import { IPipeline } from "./models";
import { normalize } from "normalizr";

export class PipelinesMutations extends Mutations<PipelinesState> {
  setActivePipelineId(id: number | null) {
    this.state.activePipelineId = id;
  }

  setEntities(payload: IPipeline[]) {
    const normalizedData = normalize<IPipeline>(payload, pipelineListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.pipelines ? normalizedData.entities.pipelines : {};
  }

  setEntity(payload: IPipeline) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  addEntity(payload: IPipeline) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = { ...this.state.entities, [payload.id]: payload };
  }

  updateEntity(payload: IPipeline) {
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
    const s = new PipelinesState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
