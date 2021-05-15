import { isEqual } from "lodash-es";
import { Mutations } from "vuex-smart-module";
import { projectListSchema, ProjectsState } from ".";
import { IProject, IProjectData } from "./models";
import { normalize } from "normalizr";

export class ProjectsMutations extends Mutations<ProjectsState> {
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
    if (!isEqual(tags, this.state.tags)) {
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

  reset() {
    // acquire initial state
    const s = new ProjectsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
