import { Mutations } from "vuex-smart-module";
import { SegmentationState } from ".";

export class SegmentationMutations extends Mutations<SegmentationState> {
  setSelectedAcquisitionIds(ids: number[]) {
    this.state.selectedAcquisitionIds = ids;
  }

  setSelectedTags(tags: string[]) {
    this.state.selectedTags = tags;
  }

  reset() {
    // acquire initial state
    const s = new SegmentationState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
