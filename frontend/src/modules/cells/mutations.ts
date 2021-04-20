import { Mutations } from "vuex-smart-module";
import { CellsState, resultListSchema } from ".";
import { ICell, ICentroidsData, IRawColorsData, IRawResultData, IRawScatterData, IResult } from "./models";
import { normalize } from "normalizr";

function getCellsByAcquisition(cells: { [p: string]: ICell }) {
  // Refresh cellsByAcquisition mapping
  const cellsByAcquisition = new Map<number, ICell[]>();
  Object.values(cells).forEach((c) => {
    if (!cellsByAcquisition.has(c.acquisitionId)) {
      cellsByAcquisition.set(c.acquisitionId, []);
    }
    cellsByAcquisition.get(c.acquisitionId)!.push(c);
  });
  return cellsByAcquisition;
}

export class CellsMutations extends Mutations<CellsState> {
  initializeCells(payload: ICentroidsData) {
    const cells = {};
    const cellsByAcquisition = new Map<number, ICell[]>();
    payload.cellIds.forEach((cellId, i) => {
      const acquisitionId = payload.acquisitionIds[i];
      if (!cellsByAcquisition.has(acquisitionId)) {
        cellsByAcquisition.set(acquisitionId, []);
      }
      const cell: ICell = {
        cellId: cellId,
        objectNumber: payload.objectNumbers[i],
        acquisitionId: acquisitionId,
        xy: [payload.x[i], payload.y[i]],
        color: payload.colors[i],
        defaultColor: payload.colors[i],
        mappings: {},
      };
      cells[cellId] = cell;
      cellsByAcquisition.get(acquisitionId)!.push(cell);
    });
    this.state.cells = Object.freeze(cells);
    this.state.cellsByAcquisition = Object.freeze(cellsByAcquisition);
  }

  updateCellsByResult(payload: IRawResultData) {
    if (!this.state.cells) {
      return;
    }
    const cells = { ...this.state.cells };

    // Reset mappings for all cells
    Object.values(cells).forEach((cell) => {
      cell.mappings = {};
    });

    for (let i = 0; i < payload.cellIds.length; i++) {
      const cell = cells[payload.cellIds[i]];
      if (payload.mappings) {
        const mappings: { [key: string]: [number, number] } = {};
        if (payload.mappings.pca) {
          mappings.pca = [payload.mappings.pca.x.data[i], payload.mappings.pca.y.data[i]];
        }
        if (payload.mappings.tsne) {
          mappings.tsne = [payload.mappings.tsne.x.data[i], payload.mappings.tsne.y.data[i]];
        }
        if (payload.mappings.umap) {
          mappings.umap = [payload.mappings.umap.x.data[i], payload.mappings.umap.y.data[i]];
        }
        cell.mappings = mappings;
      }
    }

    const cellsByAcquisition = getCellsByAcquisition(cells);

    this.state.markers = Object.freeze(payload.markers);
    this.state.cells = Object.freeze(cells);
    this.state.cellsByAcquisition = Object.freeze(cellsByAcquisition);
  }

  updateCellsByColors(payload: IRawColorsData | null) {
    if (!this.state.cells) {
      return;
    }
    const cells = { ...this.state.cells };

    // Reset color for all cells
    Object.values(cells).forEach((cell) => {
      cell.color = cell.defaultColor;
    });

    if (payload !== null) {
      for (let i = 0; i < payload.cellIds.length; i++) {
        const cell = cells[payload.cellIds[i]]!;
        cell.color = payload.colors.data[i];
      }
    }

    const cellsByAcquisition = getCellsByAcquisition(cells);

    this.state.cells = Object.freeze(cells);
    this.state.cellsByAcquisition = Object.freeze(cellsByAcquisition);
  }

  updateCellsByScatterplot(payload: IRawScatterData | null) {
    if (!this.state.cells) {
      return;
    }
    const cells = { ...this.state.cells };

    // Reset scatterplot data for all cells
    Object.values(cells).forEach((cell) => {
      delete cell.mappings.scatterplot;
    });

    if (payload !== null) {
      for (let i = 0; i < payload.cellIds.length; i++) {
        const cell = cells[payload.cellIds[i]]!;
        cell.mappings.scatterplot = [payload.x.data[i], payload.y.data[i]];
      }
    }

    const cellsByAcquisition = getCellsByAcquisition(cells);

    this.state.cells = Object.freeze(cells);
    this.state.cellsByAcquisition = Object.freeze(cellsByAcquisition);
  }

  setSelectedCellIds(payload: string[]) {
    this.state.selectedCellIds = payload;
  }

  setActiveResultId(id: number | null) {
    this.state.activeResultId = id;
    if (this.state.cells && id === null) {
      const cells = { ...this.state.cells };

      // Reset mappings for all cells
      Object.values(cells).forEach((cell) => {
        cell.mappings = {};
      });

      const cellsByAcquisition = getCellsByAcquisition(cells);

      this.state.cells = Object.freeze(cells);
      this.state.cellsByAcquisition = Object.freeze(cellsByAcquisition);
    }
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

  resetResultData() {
    this.state.selectedCellIds = [];
    this.state.heatmap = null;
    this.state.markers = [];
  }

  reset() {
    // acquire initial state
    const s = new CellsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
