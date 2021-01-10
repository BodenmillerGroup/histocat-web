import { Module } from "vuex-smart-module";
import { SegmentationActions } from "./actions";
import { SegmentationGetters } from "./getters";
import { SegmentationMutations } from "./mutations";

export class SegmentationState {
  selectedAcquisitionIds: ReadonlyArray<number> = [];
  selectedTags: string[] = [];
}

export const segmentationModule = new Module({
  namespaced: true,

  state: SegmentationState,
  getters: SegmentationGetters,
  mutations: SegmentationMutations,
  actions: SegmentationActions,
});
