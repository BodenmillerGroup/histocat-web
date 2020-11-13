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
} from "./models";
import { normalize } from "normalizr";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SELECTED_CELLS } from "./events";

export class ResultsMutations extends Mutations<ResultsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_SELECTED_CELLS, (payload) => this.setSelectedCells(payload));
  }

  setActiveResultId(id: number | null) {
    this.state.activeResultId = id;
  }

  setHeatmap(heatmap: { type: string; label: string } | null) {
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

  setCells(payload: IRawResultData) {
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
      cellData.color = payload.colors ? Number(payload.colors.data[i]) : cellData.acquisitionId;
      cells.set(cellData.cellId, cellData);
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
        y: payload.y.data[i]
      });
    }
    this.state.scatterData = Object.freeze(scatterPlotData);
  }

  setBoxPlotData(data: IPlotSeries[]) {
    this.state.boxPlotData = data;
  }

  setPhenoGraphData(data: IPhenoGraphData | null) {
    this.state.phenographData = data;
  }

  reset() {
    // acquire initial state
    const s = new ResultsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
