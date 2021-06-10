import { Getters } from "vuex-smart-module";
import { CellsState } from ".";

export class CellsGetters extends Getters<CellsState> {
  get cells() {
    return this.state.cells;
  }

  get cellsList() {
    return Object.values(this.state.cells);
  }

  get selectedCellIds() {
    return this.state.selectedCellIds;
  }

  get selectedCells() {
    return this.state.cells ? this.state.selectedCellIds.map((cellId) => this.state.cells![cellId]) : [];
  }

  get results() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  get activeResultId() {
    return this.state.activeResultId;
  }

  get activeResult() {
    return this.getters.activeResultId ? this.state.entities[this.getters.activeResultId] : null;
  }

  get markers() {
    return this.state.markers;
  }

  get heatmap() {
    return this.state.heatmap;
  }
}
