import { equals } from "rambda";
import { Mutations } from "vuex-smart-module";
import { projectListSchema, ProjectsState } from ".";
import { IProject, IProjectData } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_ACTIVE_ACQUISITION_ID,
  SET_ACTIVE_WORKSPACE_NODE,
  SET_CHANNEL_STACK_IMAGE,
  SET_SELECTED_METALS,
} from "./events";
import { normalize } from "normalizr";

export class ProjectsMutations extends Mutations<ProjectsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_ACTIVE_ACQUISITION_ID, (payload) => this.setActiveAcquisitionId(payload));
    BroadcastManager.subscribe(SET_ACTIVE_WORKSPACE_NODE, (payload) => this.setActiveWorkspaceNode(payload));
    BroadcastManager.subscribe(SET_SELECTED_METALS, (payload) => this.setSelectedMetals(payload));
    BroadcastManager.subscribe(SET_CHANNEL_STACK_IMAGE, (payload) => this.setChannelStackImage(payload));
    BroadcastManager.subscribe(SET_CHANNEL_STACK_IMAGE, (payload) => this.setChannelStackImage(payload));
  }

  setActiveProjectId(id: number | null) {
    this.state.activeProjectId = id;
  }

  setActiveAcquisitionId(id: number | null) {
    this.state.activeAcquisitionId = id;
  }

  setProjectData(data: IProjectData) {
    this.state.projectData = data;
  }

  setTags(tags: string[]) {
    if (!equals(tags, this.state.tags)) {
      this.state.tags = tags;
    }
  }

  setEntities(payload: IProject[]) {
    const normalizedData = normalize<IProject>(payload, projectListSchema);
    this.state.ids = normalizedData.result;
    this.state.entities = normalizedData.entities.projects ? Object.freeze(normalizedData.entities.projects) : {};
  }

  setEntity(payload: IProject) {
    const existingId = this.state.ids.find((id) => id === payload.id);
    if (!existingId) {
      this.state.ids = this.state.ids.concat(payload.id);
    }
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  addEntity(payload: IProject) {
    this.state.ids = this.state.ids.concat(payload.id);
    this.state.entities = Object.freeze({ ...this.state.entities, [payload.id]: payload });
  }

  updateEntity(payload: IProject) {
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
    const s = new ProjectsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
