import { IPresetCreate } from "./models";
import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { PresetsState } from ".";
import { api } from "./api";
import { PresetsGetters } from "./getters";
import { PresetsMutations } from "./mutations";
import { settingsModule } from "@/modules/settings";
import { projectsModule } from "@/modules/projects";

export class PresetsActions extends Actions<PresetsState, PresetsGetters, PresetsMutations, PresetsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  projects?: Context<typeof projectsModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.projects = projectsModule.context(store);
    this.settings = settingsModule.context(store);
  }

  async getPresets(projectId: number) {
    try {
      const data = await api.getProjectPresets(projectId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createPreset(name: string) {
    try {
      const projectId = this.projects!.getters.activeProjectId;
      const presetData = {
        channelsSettings: this.settings!.getters.channelsSettings,
      };
      const payload: IPresetCreate = {
        name: name,
        project_id: projectId!,
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
        this.settings?.actions.setPreset(preset);
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
