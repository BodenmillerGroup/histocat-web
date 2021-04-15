import { Module } from "vuex-smart-module";
import { CellsActions } from "./actions";
import { CellsGetters } from "./getters";
import { CellsMutations } from "./mutations";
import { ICell, IResult } from "./models";
import { schema } from "normalizr";

export const resultSchema = new schema.Entity("results");
export const resultListSchema = [resultSchema];

export class CellsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IResult } = {};
  activeResultId: number | null = null;

  // Map cells by acquisitionId
  cells: { [cellId: string]: ICell } | null = null;
  cellsByAcquisition: Readonly<Map<number, ICell[]>> | null = null;

  selectedCellIds: string[] = [];

  heatmap: { type: string; label: string; value: string } | null = null;
  markers: Readonly<string[]> = [];
}

export const cellsModule = new Module({
  namespaced: true,

  state: CellsState,
  getters: CellsGetters,
  mutations: CellsMutations,
  actions: CellsActions,
});
