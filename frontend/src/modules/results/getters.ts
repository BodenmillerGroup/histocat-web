import { Getters } from "vuex-smart-module";
import { ResultsState } from ".";
import { ICellData, ICellPoint } from "./models";

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

  get markers() {
    return this.state.markers;
  }

  get heatmap() {
    return this.state.heatmap;
  }

  get cells() {
    return this.state.cells;
  }

  get selectedCells() {
    return this.state.selectedCells;
  }

  get scatterData() {
    return this.state.scatterData;
  }

  get scatterPlotData() {
    if (this.getters.cells && this.getters.scatterData) {
      const output = new Map<number, ICellPoint[]>();
      this.getters.cells.forEach((value, key) => {
        if (!output.has(value.acquisitionId)) {
          output.set(value.acquisitionId, []);
        }
        const acquisitionCells = output.get(value.acquisitionId)!;
        const scatterPointData = this.getters.scatterData?.get(value.cellId)!;
        acquisitionCells.push({
          cellId: value.cellId,
          acquisitionId: value.acquisitionId,
          objectNumber: value.objectNumber,
          x: scatterPointData.x,
          y: scatterPointData.y,
          color: value.color,
        });
      });
      return output;
    }
    return null;
  }

  get cellsByAcquisition() {
    const output = new Map<number, ICellData[]>();
    this.getters.cells?.forEach((value, key) => {
      if (!output.has(value.acquisitionId)) {
        output.set(value.acquisitionId, []);
      }
      const acquisitionPoints = output.get(value.acquisitionId)!;
      acquisitionPoints.push(value);
    });
    return output;
  }

  get pcaData() {
    if (this.getters.activeResult && this.getters.activeResult.output.pca) {
      const output = new Map<number, ICellPoint[]>();
      this.getters.cells?.forEach((value, key) => {
        if (!output.has(value.acquisitionId)) {
          output.set(value.acquisitionId, []);
        }
        const acquisitionCells = output.get(value.acquisitionId)!;
        acquisitionCells.push({
          cellId: value.cellId,
          acquisitionId: value.acquisitionId,
          objectNumber: value.objectNumber,
          x: value.mappings!.pca.x,
          y: value.mappings!.pca.y,
          color: value.color,
        });
      });
      return output;
    }
    return null;
  }

  get tsneData() {
    if (this.getters.activeResult && this.getters.activeResult.output.tsne) {
      const output = new Map<number, ICellPoint[]>();
      this.getters.cells?.forEach((value, key) => {
        if (!output.has(value.acquisitionId)) {
          output.set(value.acquisitionId, []);
        }
        const acquisitionCells = output.get(value.acquisitionId)!;
        acquisitionCells.push({
          cellId: value.cellId,
          acquisitionId: value.acquisitionId,
          objectNumber: value.objectNumber,
          x: value.mappings!.tsne.x,
          y: value.mappings!.tsne.y,
          color: value.color,
        });
      });
      return output;
    }
    return null;
  }

  get umapData() {
    if (this.getters.activeResult && this.getters.activeResult.output.umap) {
      const output = new Map<number, ICellPoint[]>();
      this.getters.cells?.forEach((value, key) => {
        if (!output.has(value.acquisitionId)) {
          output.set(value.acquisitionId, []);
        }
        const acquisitionCells = output.get(value.acquisitionId)!;
        acquisitionCells.push({
          cellId: value.cellId,
          acquisitionId: value.acquisitionId,
          objectNumber: value.objectNumber,
          x: value.mappings!.umap.x,
          y: value.mappings!.umap.y,
          color: value.color,
        });
      });
      return output;
    }
    return null;
  }
}
