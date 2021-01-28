import { Getters } from "vuex-smart-module";
import { SegmentationState } from ".";

export class SegmentationGetters extends Getters<SegmentationState> {
  get selectedAcquisitionIds() {
    return this.state.selectedAcquisitionIds;
  }

  get nucleiChannels() {
    return this.state.nucleiChannels;
  }

  get cytoplasmChannels() {
    return this.state.cytoplasmChannels;
  }
}
