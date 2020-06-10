import { IPresetCreate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { PresetState } from ".";
import { api } from "./api";
import { PresetGetters } from "./getters";
import { PresetMutations } from "./mutations";

export class PresetActions extends Actions<PresetState, PresetGetters, PresetMutations, PresetActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
  }

  async getPresets(experimentId: number) {
    try {
      const data = await api.getExperimentPresets(experimentId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createPreset(payload: IPresetCreate) {
    try {
      const data = await api.createPreset(payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Preset successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getPreset(id: number) {
    try {
      const data = await api.getPreset(id);
      if (data) {
        this.mutations.setEntity(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deletePreset(id: number) {
    try {
      const data = await api.deletePreset(id);
      this.mutations.deleteEntity(data);
      this.main!.mutations.addNotification({ content: "Preset successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
