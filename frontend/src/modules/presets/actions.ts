import { IPresetCreate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { PresetState } from ".";
import { api } from "./api";
import { PresetGetters } from "./getters";
import { PresetMutations } from "./mutations";
import { settingsModule } from "@/modules/settings";
import { experimentModule } from "@/modules/experiment";

export class PresetActions extends Actions<PresetState, PresetGetters, PresetMutations, PresetActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  experiment?: Context<typeof experimentModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.experiment = experimentModule.context(store);
    this.settings = settingsModule.context(store);
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

  async createPreset(name: string) {
    try {
      const experimentId = this.experiment!.getters.activeExperimentId;
      const presetData = {
        colorMap: this.settings!.getters.colorMap,
        channelSettings: this.settings!.getters.channelSettings,
      };
      const payload: IPresetCreate = {
        name: name,
        experiment_id: experimentId!,
        data: presetData,
      };
      const data = await api.createPreset(payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Preset successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async applyPreset(id: number) {
    try {
      const preset = await api.getPreset(id);
      if (preset) {
        this.settings?.mutations.setPreset(preset);
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
