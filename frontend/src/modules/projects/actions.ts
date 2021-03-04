import { datasetsModule } from "@/modules/datasets";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { ProjectsState } from ".";
import { api } from "./api";
import { ProjectsGetters } from "./getters";
import { ExportFormat, IChannelUpdate, IProjectCreate, IProjectUpdate } from "./models";
import { ProjectsMutations } from "./mutations";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_CHANNEL_STACK_IMAGE,
  SET_ACTIVE_ACQUISITION_ID,
  SET_ACTIVE_WORKSPACE_NODE,
  SET_SELECTED_METALS,
} from "./events";
import { groupModule } from "@/modules/group";
import { resultsModule } from "@/modules/results";

export class ProjectsActions extends Actions<ProjectsState, ProjectsGetters, ProjectsMutations, ProjectsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  datasets?: Context<typeof datasetsModule>;
  results?: Context<typeof resultsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.results = resultsModule.context(store);
  }

  setActiveAcquisitionId(id?: number, isGlobal = true) {
    BroadcastManager.publish(SET_ACTIVE_ACQUISITION_ID, id, isGlobal);
  }

  setActiveWorkspaceNode(node?: { id: number; type: string }, isGlobal = true) {
    BroadcastManager.publish(SET_ACTIVE_WORKSPACE_NODE, node, isGlobal);
  }

  setSelectedMetals(metals: string[], isGlobal = true) {
    BroadcastManager.publish(SET_SELECTED_METALS, metals, isGlobal);
  }

  async getGroupProjects(groupId: number) {
    try {
      const data = await api.getGroupProjects(groupId);
      if (data) {
        this.mutations.setEntities(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getProjectsTags(groupId: number) {
    try {
      const data = await api.getProjectsTags(groupId);
      if (data) {
        this.mutations.setTags(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateProject(payload: { id: number; data: IProjectUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateProject(groupId, payload.id, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Project successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteProject(projectId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteProject(groupId, projectId);
      this.mutations.deleteEntity(projectId);
      this.main!.mutations.addNotification({ content: "Project successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createProject(payload: IProjectCreate) {
    try {
      const data = await api.createProject(payload);
      this.mutations.addEntity(data);
      this.main!.mutations.addNotification({ content: "Project successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async uploadSlide(payload: { id: number; data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      const groupId = this.group?.getters.activeGroupId!;
      await api.uploadSlide(
        this.main!.getters.token,
        groupId,
        payload.id,
        payload.data,
        () => {
          console.log("Upload has started.");
          this.main!.mutations.setProcessing(true);
        },
        () => {
          console.log("Upload completed successfully.");
          this.main!.mutations.setProcessing(false);
          this.main!.mutations.setProcessingProgress(0);
          this.main!.mutations.addNotification({ content: "File successfully uploaded", color: "success" });
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          this.main!.mutations.setProcessingProgress(percent);
        },
        () => {}
      );
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getProject(projectId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getProject(groupId, projectId);
      if (data) {
        this.mutations.setEntity(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getProjectData(projectId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getProjectData(groupId, projectId);
      if (data) {
        this.mutations.setProjectData(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStats(payload: { acquisitionId: number; channelName: string }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      return await api.getChannelStats(groupId, payload.acquisitionId, payload.channelName);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStackImage() {
    const params = await this.actions.prepareStackParams();
    if (params.channels.length === 0) {
      return;
    }
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.downloadChannelStackImage(groupId, params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        BroadcastManager.publish(SET_CHANNEL_STACK_IMAGE, reader.result);
      };
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getColorizedMaskImage() {
    const params = await this.actions.prepareStackParams();
    if (params.channels.length === 0 || !params.hasOwnProperty("mask")) {
      return;
    }
    params["mask"]["apply"] = true;
    params["mask"]["colorize"] = true;
    try {
      this.mutations.setColorizeMaskInProgress(true);
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.downloadChannelStackImage(groupId, params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        BroadcastManager.publish(SET_CHANNEL_STACK_IMAGE, reader.result);
      };
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    } finally {
      this.mutations.setColorizeMaskInProgress(false);
    }
  }

  async exportChannelStackImage(format: ExportFormat = "png") {
    const params = await this.actions.prepareStackParams(format);
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.downloadChannelStackImage(groupId, params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteSlide(id: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteSlide(groupId, id);
      this.main!.mutations.addNotification({ content: "Slide successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateChannel(payload: { acquisitionId: number; data: IChannelUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateChannel(groupId, payload.acquisitionId, payload.data);
      if (data) {
        this.mutations.setProjectData(data);
      }
      await this.actions.getChannelStackImage();
      this.main!.mutations.addNotification({ content: "Channel successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  /**
   * Prepare stack image call parameters
   */
  private prepareStackParams(format: "png" | "tiff" = "png") {
    const activeAcquisitionId = this.getters.activeAcquisitionId;
    const settings = this.settings!.getters.channelsSettings;
    const channels = this.getters.selectedChannels.map((channel) => {
      const channelSettings = settings[channel.name];
      const color = channelSettings ? channelSettings.color : undefined;
      const min = channelSettings && channelSettings.levels ? channelSettings.levels.min : undefined;
      const max = channelSettings && channelSettings.levels ? channelSettings.levels.max : undefined;
      return {
        name: channel.name,
        color: color,
        min: min,
        max: max,
      };
    });

    const filter = this.settings!.getters.filter;
    const scalebar = this.settings!.getters.scalebar;

    const output = {
      acquisitionId: activeAcquisitionId,
      format: format,
      filter: filter,
      scalebar: scalebar,
      channels: channels,
    };

    const activeDataset = this.datasets!.getters.activeDataset;
    if (activeDataset) {
      output["datasetId"] = activeDataset.id;
      const maskSettings = this.settings!.getters.maskSettings;
      const activeAcquisitionId = this.getters.activeAcquisitionId;
      if (activeAcquisitionId && activeDataset.meta.masks) {
        const mask = activeDataset.meta.masks[activeAcquisitionId];
        if (mask) {
          output["mask"] = {
            colorize: false,
            mode: maskSettings.mode,
            location: mask.location,
          };
          if (this.results?.getters.heatmap) {
            output["mask"]["colorsType"] = this.results.getters.heatmap.type;
            output["mask"]["colorsName"] = this.results.getters.heatmap.label;
          }
          if (this.results?.getters.activeResultId) {
            output["mask"]["resultId"] = this.results?.getters.activeResultId;
          }
          // Prepare selected cell ids visualisation
          const selectedCells = this.results?.getters.selectedCells?.filter(
            (v) => v.acquisitionId === activeAcquisitionId
          );
          if (selectedCells && selectedCells.length > 0) {
            output["mask"]["gated"] = true;
            output["mask"]["cellIds"] = selectedCells.map((item) => item.objectNumber);
          }
        }
      }
    }
    return output;
  }
}
