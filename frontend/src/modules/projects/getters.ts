import { Getters } from "vuex-smart-module";
import { ProjectsState } from ".";

export class ProjectsGetters extends Getters<ProjectsState> {
  get projects() {
    return this.state.ids.map((id) => this.state.entities[id]);
  }

  get tags() {
    return this.state.tags;
  }

  getProject(id: number) {
    return this.state.entities[id];
  }

  get activeProjectId() {
    return this.state.activeProjectId;
  }

  get projectData() {
    return this.state.projectData;
  }

  get activeAcquisition() {
    if (this.getters.projectData && this.getters.projectData.slides) {
      for (const slide of this.getters.projectData.slides) {
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
