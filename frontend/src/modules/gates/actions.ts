import { IGateCreate, IGateUpdate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GatesState } from ".";
import { api } from "./api";
import { GatesGetters } from "./getters";
import { GatesMutations } from "./mutations";
import { datasetsModule } from "@/modules/datasets";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_ACTIVE_GATE_ID } from "./events";
import { groupModule } from "@/modules/group";
import { resultsModule } from "@/modules/results";
import { ISelectedCell } from "@/modules/results/models";

export class GatesActions extends Actions<GatesState, GatesGetters, GatesMutations, GatesActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  dataset?: Context<typeof datasetsModule>;
  results?: Context<typeof resultsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.dataset = datasetsModule.context(store);
    this.results = resultsModule.context(store);
  }

  setActiveGateId(id: number | null, isGlobal = true) {
    BroadcastManager.publish(SET_ACTIVE_GATE_ID, id, isGlobal);
  }

  async getGates(datasetId: number) {
    try {
      const data = await api.getDatasetGates(datasetId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createGate(name: string) {
    try {
      const datasetId = this.dataset!.getters.activeDatasetId;
      const selectedCells = this.results!.getters.selectedCells;
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
        const data = await api.createGate(payload);
        this.mutations.addEntity(data);
        this.main!.mutations.addNotification({ content: "Gate successfully created", color: "success" });
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async applyGate(id: number) {
    try {
      const data = await api.getGate(id);
      if (data) {
        const selectedCells: ISelectedCell[] = [];
        for (let i = 0; i < data.acquisition_ids.length; i++) {
          const acquisitionId = data.acquisition_ids[i];
          const index = data.indices[i];
          const cellId = data.cell_ids[i];
          selectedCells.push(Object.freeze( { acquisitionId: acquisitionId, cellId: cellId, objectNumber: index }));
        }
        this.results?.actions.setSelectedCells(selectedCells);
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
      const data = await api.deleteGate(id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Gate successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
