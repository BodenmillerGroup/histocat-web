import { IGateCreate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GateState } from ".";
import { api } from "./api";
import { GateGetters } from "./getters";
import { GateMutations } from "./mutations";
import { datasetModule } from "@/modules/datasets";
import { selectionModule } from "@/modules/selection";
import { SelectedCell } from "@/modules/selection/models";

export class GateActions extends Actions<GateState, GateGetters, GateMutations, GateActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  dataset?: Context<typeof datasetModule>;
  selection?: Context<typeof selectionModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.dataset = datasetModule.context(store);
    this.selection = selectionModule.context(store);
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
      const selectedCells = this.selection!.getters.selectedCells;
      if (datasetId && selectedCells) {
        const acquisitionIds: number[] = [];
        const indices: number[] = [];
        const cellIds: number[] = [];
        selectedCells.forEach((arr, key) => {
          arr.forEach((selectedCell) => {
            acquisitionIds.push(selectedCell.acquisitionId);
            indices.push(selectedCell.index);
            cellIds.push(selectedCell.cellId);
          });
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
        const selectedCells = new Map<number, SelectedCell[]>();
        for (let i=0; i < data.acquisition_ids.length; i++) {
          const acquisitionId = data.acquisition_ids[i];
          const index = data.indices[i];
          const cellId = data.cell_ids[i];
          if (!selectedCells.has(acquisitionId)) {
            selectedCells.set(acquisitionId, []);
          }
          const ids = selectedCells.get(acquisitionId);
          ids!.push(Object.freeze(new SelectedCell(acquisitionId, index, cellId)));
        }
        this.selection?.mutations.setSelectedCells(selectedCells);
      }
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
