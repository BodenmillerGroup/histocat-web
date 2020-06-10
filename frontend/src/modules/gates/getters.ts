import { Getters } from "vuex-smart-module";
import { GateState } from ".";

export class GateGetters extends Getters<GateState> {
  get gates() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getGate(id: number) {
    return this.state.entities[id];
  }

  get activeGateId() {
    return this.state.activeGateId;
  }

  get activeGate() {
    return this.getters.activeGateId ? this.getters.getGate(this.getters.activeGateId) : null;
  }
}
