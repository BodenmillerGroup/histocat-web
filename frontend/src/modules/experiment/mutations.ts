import { equals } from "rambda";
import { Mutations } from "vuex-smart-module";
import { ExperimentState } from ".";
import { IExperiment, IShare } from "./models";

export class ExperimentMutations extends Mutations<ExperimentState> {
  setExperiments(experiments: IExperiment[]) {
    this.state.experiments = experiments;
  }

  setShares(shares: IShare[]) {
    this.state.shares = shares;
  }

  setTags(tags: string[]) {
    this.state.tags = tags;
  }

  setExperiment(experiment: IExperiment) {
    const items = this.state.experiments.filter(item => item.id !== experiment.id);
    items.push(experiment);
    this.state.experiments = items;
  }

  deleteExperiment(id: number) {
    this.state.experiments = this.state.experiments.filter(item => item.id !== id);
  }

  setSelectedMetals(metals: string[]) {
    this.state.selectedMetals = metals;
  }

  setSelectedAcquisitionIds(ids: number[]) {
    if (!equals(ids, this.state.selectedAcquisitionIds)) {
      this.state.selectedAcquisitionIds = ids;
    }
  }

  setActiveExperimentId(id?: number) {
    this.state.activeExperimentId = id;
  }

  setActiveAcquisitionId(id?: number) {
    this.state.activeAcquisitionId = id;
  }

  setActiveWorkspaceNode(node?: { id: number; type: string }) {
    this.state.activeWorkspaceNode = node;
    if (node) {
      if (node.type === "acquisition") {
        this.setActiveAcquisitionId(node.id);
      }
    }
  }

  setChannelStackImage(base64Image: string | ArrayBuffer | null) {
    this.state.channelStackImage = base64Image;
  }

  setColorizeMaskInProgress(status: boolean) {
    this.state.colorizeMaskInProgress = status;
  }

  reset() {
    this.state.activeWorkspaceNode = undefined;
    this.state.activeExperimentId = undefined;
    this.state.activeAcquisitionId = undefined;
    this.state.selectedAcquisitionIds = [];
    this.state.selectedMetals = [];
    this.state.channelStackImage = null;
    this.state.colorizeMaskInProgress = false;
  }
}
