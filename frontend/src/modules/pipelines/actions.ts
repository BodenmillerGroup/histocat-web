import { IPipelineCreate, IPipelineUpdate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { PipelinesState } from ".";
import { api } from "./api";
import { PipelinesGetters } from "./getters";
import { PipelinesMutations } from "./mutations";
import { settingsModule } from "@/modules/settings";
import { projectsModule } from "@/modules/projects";
import { groupModule } from "@/modules/group";
import { datasetsModule } from "@/modules/datasets";

export class PipelinesActions extends Actions<PipelinesState, PipelinesGetters, PipelinesMutations, PipelinesActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  projects?: Context<typeof projectsModule>;
  datasets?: Context<typeof datasetsModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.projects = projectsModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.settings = settingsModule.context(store);
  }

  async getPipelines(projectId: number) {
    try {
      const data = await api.getProjectPipelines(projectId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createPipeline(name: string) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const projectId = this.projects!.getters.activeProjectId;
      const steps = this.getters.steps;
      const payload: IPipelineCreate = {
        name: name,
        project_id: projectId!,
        steps: steps,
      };
      const data = await api.createPipeline(groupId, payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Pipeline successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updatePipeline(payload: { pipelineId: number; data: IPipelineUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updatePipeline(groupId, payload.pipelineId, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Pipeline successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async loadPipeline(id: number) {
    try {
      const pipeline = await api.getPipeline(id);
      if (pipeline) {
        this.mutations.setSteps(pipeline.steps);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deletePipeline(id: number) {
    try {
      const data = await api.deletePipeline(id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Pipeline successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async processPipeline() {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const datasetId = this.datasets?.getters.activeDatasetId!;
      const selectedAcquisitionsIds = this.getters.selectedAcquisitionIds;
      const steps = this.getters.steps;
      const data = await api.processPipeline(groupId, {
        dataset_id: datasetId,
        acquisition_ids: selectedAcquisitionsIds,
        steps: steps,
      });
      this.main!.mutations.addNotification({ content: "Pipeline successfully submitted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
