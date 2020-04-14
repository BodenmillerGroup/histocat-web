import { Getters } from "vuex-smart-module";
import { ExperimentState } from ".";

export class ExperimentGetters extends Getters<ExperimentState> {
  get experiments() {
    return this.state.experiments;
  }

  get shares() {
    return this.state.shares;
  }

  getExperiment(id?: number) {
    return this.getters.experiments.find((item) => item.id === id);
  }

  get activeExperiment() {
    return this.getters.getExperiment(this.getters.activeExperimentId);
  }

  get activeAcquisition() {
    if (this.getters.activeExperiment && this.getters.activeExperiment.slides) {
      for (const slide of this.getters.activeExperiment.slides) {
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

  get tags() {
    return this.state.tags;
  }

  get selectedAcquisitionIds() {
    return this.state.selectedAcquisitionIds;
  }

  get activeExperimentId() {
    return this.state.activeExperimentId;
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
