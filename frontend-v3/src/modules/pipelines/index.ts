import create from "zustand";
import { schema } from "normalizr";
import { IPipeline } from "./models";

export const datasetSchema = new schema.Entity("datasets");
export const datasetListSchema = [datasetSchema];

type PipelinesState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IPipeline };
  activePipelineId: number | null;
  selectedAcquisitionIds: number[];
  steps: any[];

  setActivePipelineId(id: number | null): void;
  getActivePipeline(): IPipeline | null;
  setSteps(steps: any[]): void;
  setSelectedAcquisitionIds(ids: number[]): void;
};

export const usePipelinesStore = create<PipelinesState>((set, get) => ({
  ids: [],
  entities: {},
  activePipelineId: null,
  selectedAcquisitionIds: [],
  steps: [],

  setActivePipelineId(id: number | null) {
    set({ activePipelineId: id });
  },

  getActivePipeline() {
    const activePipelineId = get().activePipelineId;
    return activePipelineId ? get().entities[activePipelineId] : null;
  },

  setSteps(steps: any[]) {
    set({ steps: steps });
  },

  setSelectedAcquisitionIds(ids: number[]) {
    set({ selectedAcquisitionIds: ids });
  },
}));