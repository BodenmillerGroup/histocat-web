import create from "zustand";
import { schema, normalize } from "normalizr";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { AppToaster } from "../../utils/toaster";
import { useGroupsStore } from "../groups";
import {
  ICellData,
  IRawColorsData,
  IRawResultData,
  IRawScatterData,
  IResult,
  IResultUpdate,
  ISelectedCell,
} from "./models";
import { usePipelinesStore } from "../pipelines";

export const resultSchema = new schema.Entity("results");
export const resultListSchema = [resultSchema];

type ResultsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IResult };
  activeResultId: number | null;
  heatmap: { type: string; label: string } | null;
  markers: Readonly<string[]>;
  cells: Readonly<Map<string, ICellData>> | null;
  selectedCells: ISelectedCell[];
  scatterData: Readonly<Map<string, { x: number; y: number }>> | null;

  setActiveResultId(id: number | null): void;
  resetResultData(): void;
  setData(payload: IRawResultData): void;
  getActiveResult(): IResult | null;
  setSelectedCells(cells: ISelectedCell[], isGlobal?: boolean): void;
  getDatasetResults(datasetId: number): Promise<void>;
  getResultData(resultId: number): Promise<void>;
  setColors(payload: IRawColorsData | null): void;
  getColorsData(): Promise<void>;
  setScatterData(params: IRawScatterData | null): void;
  getScatterPlotData(markerX: string, markerY: string): Promise<void>;
  updateResult(resultId: number, params: IResultUpdate): Promise<void>;
  deleteResult(resultId: number): Promise<void>;
};

export const useResultsStore = create<ResultsState>((set, get) => ({
  ids: [],
  entities: {},
  activeResultId: null,
  heatmap: null,
  markers: [],
  cells: null,
  selectedCells: [],
  scatterData: null,

  setActiveResultId(id: number | null) {
    set({ activeResultId: id });
  },

  resetResultData() {
    set({
      heatmap: null,
      markers: [],
      cells: null,
      selectedCells: [],
      scatterData: null,
    });
  },

  setData(params: IRawResultData) {
    const cells = new Map<string, ICellData>();
    for (let i = 0; i < params.cellIds.length; i++) {
      const cellData: ICellData = {
        cellId: params.cellIds[i],
        acquisitionId: params.acquisitionIds[i],
        objectNumber: params.objectNumbers[i],
        x: params.x[i],
        y: params.y[i],
      };
      if (params.mappings) {
        const mappings: { [key: string]: { x: number; y: number } } = {};
        if (params.mappings.pca) {
          mappings.pca = { x: params.mappings.pca.x.data[i], y: params.mappings.pca.y.data[i] };
        }
        if (params.mappings.tsne) {
          mappings.tsne = { x: params.mappings.tsne.x.data[i], y: params.mappings.tsne.y.data[i] };
        }
        if (params.mappings.umap) {
          mappings.umap = { x: params.mappings.umap.x.data[i], y: params.mappings.umap.y.data[i] };
        }
        cellData.mappings = mappings;
      }
      cellData.color = cellData.acquisitionId;
      cells.set(cellData.cellId, cellData);
    }
    set({ markers: Object.freeze(params.markers), cells: Object.freeze(cells) });
  },

  getActiveResult() {
    const activeResultId = get().activeResultId;
    return activeResultId ? get().entities[activeResultId] : null;
  },

  setSelectedCells(cells: ISelectedCell[], isGlobal = true) {
    set({ selectedCells: cells });
  },

  async getDatasetResults(datasetId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getDatasetResults(groupId, datasetId);
      if (data) {
        const normalizedData = normalize<IResult>(data, resultListSchema);
        set({
          ids: normalizedData.result,
          entities: normalizedData.entities.results ? normalizedData.entities.results : {},
        });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getResultData(resultId: number) {
    try {
      get().setActiveResultId(resultId);
      get().resetResultData();
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getResultData(groupId, resultId);
      get().setData(data);
      const result = get().getActiveResult();
      if (result) {
        usePipelinesStore.getState().setSteps(result.pipeline);
        usePipelinesStore.getState().setSelectedAcquisitionIds(result.input);
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  setColors(params: IRawColorsData | null) {
    if (!get().cells) {
      return;
    }
    const cells = new Map<string, ICellData>();
    if (params === null) {
      get().cells!.forEach((v, k) => {
        v.color = v.acquisitionId;
        cells.set(k, v);
      });
    } else {
      for (let i = 0; i < params.cellIds.length; i++) {
        const cellData = get().cells!.get(params.cellIds[i])!;
        cellData.color = Number(params.colors.data[i]);
        cells.set(cellData.cellId, cellData);
      }
    }
    set({ cells: Object.freeze(cells) });
  },

  async getColorsData() {
    try {
      const heatmap = get().heatmap;
      const colorsType = heatmap ? heatmap.type : undefined;
      const colorsName = heatmap ? heatmap.label : undefined;
      if (!colorsType || !colorsName) {
        get().setColors(null);
        return;
      }
      const groupId = useGroupsStore.getState().activeGroupId!;
      const resultId = get().activeResultId!;
      const data = await api.getColorsData(groupId, resultId, colorsType, colorsName);
      get().setColors(data);
    } catch (error) {
      displayApiError(error);
    }
  },

  setScatterData(params: IRawScatterData | null) {
    if (params == null) {
      set({ scatterData: null });
      return;
    }
    const scatterPlotData = new Map<string, { x: number; y: number }>();
    for (let i = 0; i < params.cellIds.length; i++) {
      scatterPlotData.set(params.cellIds[i], {
        x: params.x.data[i],
        y: params.y.data[i],
      });
    }
    set({ scatterData: Object.freeze(scatterPlotData) });
  },

  async getScatterPlotData(markerX: string, markerY: string) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const resultId = get().activeResultId!;
      const data = await api.getScatterPlotData(groupId, resultId, markerX, markerY);
      if (data) {
        get().setScatterData(data);
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateResult(resultId: number, params: IResultUpdate) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.updateResult(groupId, resultId, params);
      if (data) {
        set({ entities: { ...get().entities, [data.id]: data } });
        AppToaster.show({ message: "Result successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteResult(resultId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.deleteResult(groupId, resultId);
      if (data) {
        const entities = { ...get().entities };
        delete entities[resultId];
        set({ ids: get().ids.filter((item) => item !== resultId), entities: entities });
        AppToaster.show({ message: "Result successfully deleted", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },
}));
