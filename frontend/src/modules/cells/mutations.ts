import { Mutations } from "vuex-smart-module";
import { CellsState, resultListSchema } from ".";
import { ICell, ICentroidsData, IRawColorsData, IRawResultData, IRawScatterData, IResult } from "./models";
import { normalize } from "normalizr";
import { IAnnotation } from "@/modules/annotations/models";

export class CellsMutations extends Mutations<CellsState> {
  initializeCells(payload: ICentroidsData) {
    const cells = {};
    payload.cellIds.forEach((cellId, i) => {
      const acquisitionId = payload.acquisitionIds[i];
      const cell: ICell = {
        index: i,
        cellId: cellId,
        objectNumber: payload.objectNumbers[i],
        acquisitionId: acquisitionId,
        xy: [payload.x[i], payload.y[i]],
        color: payload.colors[i],
        defaultColor: payload.colors[i],
        mappings: {},
      };
      cells[cellId] = cell;
    });
    this.state.cells = Object.freeze(cells);
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

    this.state.markers = Object.freeze(payload.markers);
    this.state.cells = Object.freeze(cells);
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

    this.state.cells = Object.freeze(cells);
  }

  updateCellsByAnnotations(payload: { annotations: IAnnotation[]; cellClasses: { [name: string]: string } }) {
    if (!this.state.cells) {
      return;
    }
    const cells = { ...this.state.cells };

    // Reset color for all cells
    Object.values(cells).forEach((cell) => {
      cell.color = cell.defaultColor;
    });

    let cellIds: string[] = [];
    let colors: string[] = [];
    payload.annotations.forEach((annotation) => {
      if (annotation.visible) {
        cellIds = cellIds.concat(annotation.cellIds);
        colors = colors.concat(Array(cellIds.length).fill(payload.cellClasses[annotation.cellClass]));
      }
    });

    for (let i = 0; i < cellIds.length; i++) {
      const cell = cells[cellIds[i]];
      cell.color = colors[i];
    }

    this.state.cells = Object.freeze(cells);
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

    this.state.cells = Object.freeze(cells);
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

      this.state.cells = Object.freeze(cells);
    }
  }

  setHeatmap(heatmap: { type: string; label: string; value: string } | null) {
    this.state.heatmap = heatmap;
  }

  setMarkers(markers: string[]) {
    this.state.markers = markers;
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
