import { Getters } from "vuex-smart-module";
import { PipelinesState } from ".";

export class PipelinesGetters extends Getters<PipelinesState> {
  get pipelines() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getPipeline(id: number) {
    return this.state.entities[id];
  }

  get activePipelineId() {
    return this.state.activePipelineId;
  }

  get activePipeline() {
    return this.getters.activePipelineId ? this.getters.getPipeline(this.getters.activePipelineId) : null;
  }
}
