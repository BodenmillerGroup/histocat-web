import { Getters } from 'vuex-smart-module';
import { ExperimentState } from '.';

export class ExperimentGetters extends Getters<ExperimentState> {
  get experiments() {
    return this.state.experiments;
  }

  get datasets() {
    return this.state.datasets;
  }

  get activeExperiment() {
    return this.state.experiments.find((item) => item.id === this.state.activeExperimentId);
  }

  get activeAcquisition() {
    if (this.activeExperiment && this.activeExperiment.slides) {
      for (const slide of this.activeExperiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            const acquisition = roi.acquisitions.find((item) => item.id === this.state.activeAcquisitionId);
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
      return this.activeAcquisition.channels.filter((channel) => {
        if (this.state.selectedMetals.includes(channel.metal)) {
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

  get activeSlideId() {
    return this.state.activeSlideId;
  }

  get activePanoramaId() {
    return this.state.activePanoramaId;
  }

  adminOneExperiment(id: number) {
    return this.state.experiments.find((item) => item.id === id);
  }

  get activeWorkspaceNode() {
    return this.state.activeWorkspaceNode;
  }

  get channelStackImage() {
    return this.state.channelStackImage;
  }
}
