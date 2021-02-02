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
import { groupModule } from "@/modules/group";

export class PresetsActions extends Actions<PresetsState, PresetsGetters, PresetsMutations, PresetsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  projects?: Context<typeof projectsModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
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
      const groupId = this.group?.getters.activeGroupId!;
      const projectId = this.projects!.getters.activeProjectId;
      const presetData = {
        channelsSettings: this.settings!.getters.channelsSettings,
        selectedTags: this.projects!.getters.selectedMetals,
      };
      const payload: IPresetCreate = {
        name: name,
        project_id: projectId!,
        data: presetData,
      };
      const data = await api.createPreset(groupId, payload);
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
        if (preset.data["channelsSettings"]) {
          this.settings?.actions.setChannelSettings(preset.data["channelsSettings"]);
        }
        if (preset.data["selectedTags"]) {
          this.projects?.actions.setSelectedMetals(preset.data["selectedTags"]);
        }
        this.projects!.actions.getChannelStackImage();
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
