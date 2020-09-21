import { Getters } from "vuex-smart-module";
import { PresetsState } from ".";

export class PresetsGetters extends Getters<PresetsState> {
  get presets() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getPreset(id: number) {
    return this.state.entities[id];
  }

  get activePresetId() {
    return this.state.activePresetId;
  }

  get activePreset() {
    return this.getters.activePresetId ? this.getters.getPreset(this.getters.activePresetId) : null;
  }
}
