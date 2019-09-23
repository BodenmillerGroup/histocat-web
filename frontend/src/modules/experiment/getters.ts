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
    return this.experiments.find(item => item.id === id);
  }

  get activeExperiment() {
    return this.getExperiment(this.activeExperimentId);
  }

  get activeAcquisition() {
    if (this.activeExperiment && this.activeExperiment.slides) {
      for (const slide of this.activeExperiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            const acquisition = roi.acquisitions.find(item => item.id === this.state.activeAcquisitionId);
            if (acquisition) {
              return acquisition;
            }
          }
        }
      }
    }
    return undefined;
  }

  get selectedMetals() {
    return this.state.selectedMetals;
  }

  get selectedChannels() {
    if (this.activeAcquisition) {
      return this.activeAcquisition.channels.filter(channel => {
        if (this.selectedMetals.includes(channel.metal)) {
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
