import { Getters } from "vuex-smart-module";
import { ExperimentState } from ".";

export class ExperimentGetters extends Getters<ExperimentState> {
  get experiments() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  get tags() {
    return this.state.tags;
  }

  getExperiment(id: number) {
    return this.state.entities[id];
  }

  get activeExperimentId() {
    return this.state.activeExperimentId;
  }

  get activeExperiment() {
    return this.getters.activeExperimentId ? this.getters.getExperiment(this.getters.activeExperimentId) : null;
  }

  get experimentData() {
    return this.state.experimentData;
  }

  get activeAcquisition() {
    if (this.getters.experimentData && this.getters.experimentData.slides) {
      for (const slide of this.getters.experimentData.slides) {
        const acquisition = slide.acquisitions.find((item) => item.id === this.state.activeAcquisitionId);
        if (acquisition) {
          return acquisition;
        }
      }
    }
    return undefined;
  }

  get selectedMetals() {
    return this.state.selectedMetals;
  }

  get selectedChannels() {
    if (this.getters.activeAcquisition) {
      return Object.values(this.getters.activeAcquisition.channels).filter((channel) => {
        if (this.getters.selectedMetals.includes(channel.name)) {
          return channel;
        }
      });
    } else {
      return [];
    }
  }

  get selectedAcquisitionIds() {
    return this.state.selectedAcquisitionIds;
  }

  get activeAcquisitionId() {
    return this.state.activeAcquisitionId;
  }

  get activeWorkspaceNode() {
    return this.state.activeWorkspaceNode;
  }

  get channelStackImage() {
    return this.state.channelStackImage;
  }

  get colorizeMaskInProgress() {
    return this.state.colorizeMaskInProgress;
  }
}
