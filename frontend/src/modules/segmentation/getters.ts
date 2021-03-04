import { Getters } from "vuex-smart-module";
import { SegmentationState } from ".";

export class SegmentationGetters extends Getters<SegmentationState> {
  get selectedAcquisitionIds() {
    return this.state.selectedAcquisitionIds;
  }

  get channels() {
    return this.state.channels;
  }

  get nucleiChannels() {
    return this.state.nucleiChannels;
  }

  get cytoplasmChannels() {
    return this.state.cytoplasmChannels;
  }
}
