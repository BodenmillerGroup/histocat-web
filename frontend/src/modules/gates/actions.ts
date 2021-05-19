import { IGateCreate, IGateUpdate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GatesState } from ".";
import { api } from "./api";
import { GatesGetters } from "./getters";
import { GatesMutations } from "./mutations";
import { datasetsModule } from "@/modules/datasets";
import { groupModule } from "@/modules/group";
import { cellsModule } from "@/modules/cells";
import { annotationsModule } from "@/modules/annotations";
import { projectsModule } from "@/modules/projects";
import { uiModule } from "@/modules/ui";

export class GatesActions extends Actions<GatesState, GatesGetters, GatesMutations, GatesActions> {
  // Declare context type
  ui?: Context<typeof uiModule>;
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  datasets?: Context<typeof datasetsModule>;
  cells?: Context<typeof cellsModule>;
  annotations?: Context<typeof annotationsModule>;
  projects?: Context<typeof projectsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.ui = uiModule.context(store);
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.cells = cellsModule.context(store);
    this.annotations = annotationsModule.context(store);
    this.projects = projectsModule.context(store);
  }

  async getGates() {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const datasetId = this.datasets!.getters.activeDatasetId;
      if (datasetId) {
        const data = await api.getDatasetGates(groupId, datasetId);
        if (data) {
          this.mutations.setEntities(data);
        }
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createGate(name: string) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const datasetId = this.datasets!.getters.activeDatasetId;
      const cellClasses = this.annotations?.getters.cellClasses;
      const annotations = this.annotations?.getters.annotations;
      if (datasetId && cellClasses && annotations) {
        const payload: IGateCreate = {
          dataset_id: datasetId!,
          name: name,
          cell_classes: cellClasses,
          annotations: annotations,
        };
        const data = await api.createGate(groupId, payload);
        this.mutations.addEntity(data);
        this.main!.mutations.addNotification({ content: "Gate successfully created", color: "success" });
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async loadGate(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getGate(groupId, id);
      if (data) {
        this.annotations?.mutations.setCellClasses(data.cell_classes);
        this.annotations?.mutations.setAnnotations(data.annotations);

        if (this.cells!.getters.heatmap && this.cells!.getters.heatmap.type === "annotation") {
          this.projects!.actions.getAnnotationData();
          if (this.ui!.getters.showMask) {
            await this.projects!.actions.getChannelStackImage();
          }
        }
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateGate(payload: { gateId: number; data: IGateUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateGate(groupId, payload.gateId, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Gate successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteGate(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteGate(groupId, id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Gate successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
