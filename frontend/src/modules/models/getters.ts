import { Getters } from "vuex-smart-module";
import { ModelsState } from ".";

export class ModelsGetters extends Getters<ModelsState> {
  get models() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getModel(id: number) {
    return this.state.entities[id];
  }

  get activeModelId() {
    return this.state.activeModelId;
  }

  get activeModel() {
    return this.getters.activeModelId ? this.getters.getModel(this.getters.activeModelId) : null;
  }
}
