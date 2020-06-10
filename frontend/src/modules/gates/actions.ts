import { IGateCreate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GateState } from ".";
import { api } from "./api";
import { GateGetters } from "./getters";
import { GateMutations } from "./mutations";

export class GateActions extends Actions<GateState, GateGetters, GateMutations, GateActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
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

  async createGate(payload: IGateCreate) {
    try {
      const data = await api.createGate(payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Gate successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getGate(id: number) {
    try {
      const data = await api.getGate(id);
      if (data) {
        this.mutations.setEntity(data);
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
