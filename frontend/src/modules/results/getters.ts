import { Getters } from "vuex-smart-module";
import { ResultsState } from ".";

export class ResultsGetters extends Getters<ResultsState> {
  get results() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  getResult(id: number) {
    return this.state.entities[id];
  }

  get activeResultId() {
    return this.state.activeResultId;
  }

  get activeResult() {
    return this.getters.activeResultId ? this.getters.getResult(this.getters.activeResultId) : null;
  }

  get heatmap() {
    return this.state.heatmap;
  }
}
