import { Mutations } from "vuex-smart-module";
import { SegmentationState } from ".";

export class SegmentationMutations extends Mutations<SegmentationState> {
  setSelectedAcquisitionIds(ids: number[]) {
    this.state.selectedAcquisitionIds = ids;
  }

  setNucleiChannels(channels: string[]) {
    this.state.nucleiChannels = channels;
  }

  setCytoplasmChannels(channels: string[]) {
    this.state.cytoplasmChannels = channels;
  }

  reset() {
    // acquire initial state
    const s = new SegmentationState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
