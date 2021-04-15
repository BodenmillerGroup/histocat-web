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

export class GatesActions extends Actions<GatesState, GatesGetters, GatesMutations, GatesActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  dataset?: Context<typeof datasetsModule>;
  cells?: Context<typeof cellsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.dataset = datasetsModule.context(store);
    this.cells = cellsModule.context(store);
  }

  async getGates(datasetId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getDatasetGates(groupId, datasetId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createGate(name: string) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const datasetId = this.dataset!.getters.activeDatasetId;
      const selectedCells = this.cells!.getters.selectedCells;
      if (datasetId && selectedCells) {
        const acquisitionIds: number[] = [];
        const indices: number[] = [];
        const cellIds: string[] = [];
        selectedCells.forEach((selectedCell) => {
          acquisitionIds.push(selectedCell.acquisitionId);
          indices.push(selectedCell.objectNumber);
          cellIds.push(selectedCell.cellId);
        });
        const payload: IGateCreate = {
          dataset_id: datasetId!,
          name: name,
          acquisition_ids: acquisitionIds,
          indices: indices,
          cell_ids: cellIds,
        };
        const data = await api.createGate(groupId, payload);
        this.mutations.addEntity(data);
        this.main!.mutations.addNotification({ content: "Gate successfully created", color: "success" });
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async applyGate(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getGate(groupId, id);
      if (data) {
        const selectedCellIds: string[] = [];
        for (let i = 0; i < data.acquisition_ids.length; i++) {
          // const acquisitionId = data.acquisition_ids[i];
          // const index = data.indices[i];
          const cellId = data.cell_ids[i];
          selectedCellIds.push(cellId);
        }
        this.cells?.mutations.setSelectedCellIds(selectedCellIds);
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
