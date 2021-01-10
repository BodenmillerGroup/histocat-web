import { Getters } from "vuex-smart-module";
import { SegmentationState } from ".";

export class SegmentationGetters extends Getters<SegmentationState> {
  get selectedAcquisitionIds() {
    return this.state.selectedAcquisitionIds;
  }

  selectedTags(id: number) {
    return this.state.selectedTags;
  }
}
