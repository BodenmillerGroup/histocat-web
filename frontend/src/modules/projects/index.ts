import { Module } from "vuex-smart-module";
import { ProjectsActions } from "./actions";
import { ProjectsGetters } from "./getters";
import { IProject, IProjectData } from "./models";
import { ProjectsMutations } from "./mutations";
import { schema } from "normalizr";

export const projectSchema = new schema.Entity("projects");
export const projectListSchema = [projectSchema];

export class ProjectsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IProject } = {};
  tags: string[] = [];
  projectData: IProjectData | null = null;

  activeProjectId: number | null = null;
  activeAcquisitionId: number | null = null;
  activeWorkspaceNode: object | null = null;
  selectedAcquisitionIds: number[] = [];
  selectedMetals: string[] = [];
  channelStackImage: string | ArrayBuffer | null = null;

  colorizeMaskInProgress = false;
}

export const projectsModule = new Module({
  namespaced: true,

  state: ProjectsState,
  getters: ProjectsGetters,
  mutations: ProjectsMutations,
  actions: ProjectsActions,
});
