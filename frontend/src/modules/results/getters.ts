import { Getters } from "vuex-smart-module";
import { ResultState } from ".";

export class ResultGetters extends Getters<ResultState> {
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

  get tsneResults() {
    return this.getters.results.filter((v) => v.type === "tsne");
  }

  get umapResults() {
    return this.getters.results.filter((v) => v.type === "umap");
  }

  get phenographResults() {
    return this.getters.results.filter((v) => v.type === "phenograph");
  }
}
