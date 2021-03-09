import { Mutations } from "vuex-smart-module";
import { resultListSchema, ResultsState } from ".";
import {
  ICellData,
  IPhenoGraphData,
  IPlotSeries,
  IRawResultData,
  IResult,
  IRawScatterData,
  ISelectedCell,
  IRawColorsData,
} from "./models";
import { normalize } from "normalizr";

export class ResultsMutations extends Mutations<ResultsState> {
  setActiveResultId(id: number | null) {
    this.state.activeResultId = id;
  }

  setHeatmap(heatmap: { type: string; label: string; value: string } | null) {
    this.state.heatmap = heatmap;
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

  setData(payload: IRawResultData) {
    const cells = new Map<string, ICellData>();
    for (let i = 0; i < payload.cellIds.length; i++) {
      const cellData: ICellData = {
        cellId: payload.cellIds[i],
        acquisitionId: payload.acquisitionIds[i],
        objectNumber: payload.objectNumbers[i],
        x: payload.x[i],
        y: payload.y[i],
      };
      if (payload.mappings) {
        const mappings: { [key: string]: { x: number; y: number } } = {};
        if (payload.mappings.pca) {
          mappings.pca = { x: payload.mappings.pca.x.data[i], y: payload.mappings.pca.y.data[i] };
        }
        if (payload.mappings.tsne) {
          mappings.tsne = { x: payload.mappings.tsne.x.data[i], y: payload.mappings.tsne.y.data[i] };
        }
        if (payload.mappings.umap) {
          mappings.umap = { x: payload.mappings.umap.x.data[i], y: payload.mappings.umap.y.data[i] };
        }
        cellData.mappings = mappings;
      }
      cellData.color = cellData.acquisitionId;
      cells.set(cellData.cellId, cellData);
    }
    this.state.markers = Object.freeze(payload.markers);
    this.state.cells = Object.freeze(cells);
  }

  setColors(payload: IRawColorsData | null) {
    if (!this.state.cells) {
      return;
    }
    const cells = new Map<string, ICellData>();
    if (payload === null) {
      this.state.cells.forEach((v, k) => {
        v.color = v.acquisitionId;
        cells.set(k, v);
      });
    } else {
      for (let i = 0; i < payload.cellIds.length; i++) {
        const cellData = this.state.cells.get(payload.cellIds[i])!;
        cellData.color = Number(payload.colors.data[i]);
        cells.set(cellData.cellId, cellData);
      }
    }
    this.state.cells = Object.freeze(cells);
  }

  setSelectedCells(payload: ISelectedCell[]) {
    this.state.selectedCells = payload;
  }

  setScatterData(payload: IRawScatterData | null) {
    if (payload == null) {
      this.state.scatterData = null;
      return;
    }
    const scatterPlotData = new Map<string, { x: number; y: number }>();
    for (let i = 0; i < payload.cellIds.length; i++) {
      scatterPlotData.set(payload.cellIds[i], {
        x: payload.x.data[i],
        y: payload.y.data[i],
      });
    }
    this.state.scatterData = Object.freeze(scatterPlotData);
  }

  resetResultData() {
    this.state.heatmap = null;
    this.state.markers = [];
    this.state.cells = null;
    this.state.selectedCells = [];
    this.state.scatterData = null;
  }

  reset() {
    // acquire initial state
    const s = new ResultsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
