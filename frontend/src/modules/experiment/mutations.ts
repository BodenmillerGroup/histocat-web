import { equals } from "rambda";
import { Mutations } from "vuex-smart-module";
import { ExperimentState } from ".";
import { IExperiment, IShare } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_ACTIVE_ACQUISITION_ID,
  SET_ACTIVE_WORKSPACE_NODE,
  SET_CHANNEL_STACK_IMAGE,
  SET_SELECTED_METALS,
} from "@/modules/experiment/events";

export class ExperimentMutations extends Mutations<ExperimentState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_ACTIVE_ACQUISITION_ID, (payload) => this.setActiveAcquisitionId(payload));
    BroadcastManager.subscribe(SET_ACTIVE_WORKSPACE_NODE, (payload) => this.setActiveWorkspaceNode(payload));
    BroadcastManager.subscribe(SET_SELECTED_METALS, (payload) => this.setSelectedMetals(payload));
    BroadcastManager.subscribe(SET_CHANNEL_STACK_IMAGE, (payload) => this.setChannelStackImage(payload));
  }

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
    const items = this.state.experiments.filter((item) => item.id !== experiment.id);
    items.push(experiment);
    this.state.experiments = items;
  }

  deleteExperiment(id: number) {
    this.state.experiments = this.state.experiments.filter((item) => item.id !== id);
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
  }

  setChannelStackImage(base64Image: string | ArrayBuffer | null) {
    this.state.channelStackImage = base64Image;
  }

  setColorizeMaskInProgress(status: boolean) {
    this.state.colorizeMaskInProgress = status;
  }

  reset() {
    // acquire initial state
    const s = new ExperimentState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
