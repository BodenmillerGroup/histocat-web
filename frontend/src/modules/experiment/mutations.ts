import { equals } from "rambda";
import { Mutations } from "vuex-smart-module";
import { experimentListSchema, ExperimentState } from ".";
import { IExperiment, IExperimentData } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_ACTIVE_ACQUISITION_ID,
  SET_ACTIVE_WORKSPACE_NODE,
  SET_CHANNEL_STACK_IMAGE,
  SET_SELECTED_ACQUISITION_IDS,
  SET_SELECTED_METALS,
} from "./events";
import { normalize } from "normalizr";

export class ExperimentMutations extends Mutations<ExperimentState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_ACTIVE_ACQUISITION_ID, (payload) => this.setActiveAcquisitionId(payload));
    BroadcastManager.subscribe(SET_ACTIVE_WORKSPACE_NODE, (payload) => this.setActiveWorkspaceNode(payload));
    BroadcastManager.subscribe(SET_SELECTED_METALS, (payload) => this.setSelectedMetals(payload));
    BroadcastManager.subscribe(SET_CHANNEL_STACK_IMAGE, (payload) => this.setChannelStackImage(payload));
    BroadcastManager.subscribe(SET_CHANNEL_STACK_IMAGE, (payload) => this.setChannelStackImage(payload));
    BroadcastManager.subscribe(SET_SELECTED_ACQUISITION_IDS, (payload) => this.setSelectedAcquisitionIds(payload));
  }

  setActiveExperimentId(id: number | null) {
    this.state.activeExperimentId = id;
  }

  setActiveAcquisitionId(id: number | null) {
    this.state.activeAcquisitionId = id;
  }

  setExperimentData(data: IExperimentData) {
    this.state.experimentData = data;
  }

  setTags(tags: string[]) {
    if (!equals(tags, this.state.tags)) {
      this.state.tags = tags;
    }
  }

  setEntities(payload: IExperiment[]) {
    const normalizedData = normalize<IExperiment>(payload, experimentListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.experiments ? Object.freeze(normalizedData.entities.experiments) : {};
  }

  setEntity(payload: IExperiment) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  addEntity(payload: IExperiment) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  updateEntity(payload: IExperiment) {
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  deleteEntity(id: number) {
    this.state.ids = this.state.ids.filter((item) => item !== id);
    const entities = { ...this.state.entities };
    delete entities[id];
    this.state.entities = Object.freeze(entities);
  }

  setSelectedMetals(metals: string[]) {
    this.state.selectedMetals = metals;
  }

  setSelectedAcquisitionIds(ids: number[]) {
    if (!equals(ids, this.state.selectedAcquisitionIds)) {
      this.state.selectedAcquisitionIds = ids;
    }
  }

  setActiveWorkspaceNode(node: { id: number; type: string } | null) {
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
