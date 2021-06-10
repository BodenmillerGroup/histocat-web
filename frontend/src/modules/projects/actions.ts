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
import { groupModule } from "@/modules/group";
import { cellsModule } from "@/modules/cells";
import { annotationsModule } from "@/modules/annotations";
import { uiModule } from "@/modules/ui";

export class ProjectsActions extends Actions<ProjectsState, ProjectsGetters, ProjectsMutations, ProjectsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  ui?: Context<typeof uiModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  datasets?: Context<typeof datasetsModule>;
  cells?: Context<typeof cellsModule>;
  annotations?: Context<typeof annotationsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.ui = uiModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.cells = cellsModule.context(store);
    this.annotations = annotationsModule.context(store);
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
          this.ui!.mutations.setProcessing(true);
        },
        () => {
          console.log("Upload completed successfully.");
          this.ui!.mutations.setProcessing(false);
          this.ui!.mutations.setProcessingProgress(0);
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          this.ui!.mutations.setProcessingProgress(percent);
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
      const data = await api.downloadChannelStackImage(groupId, params);
      if (data) {
        const blob = await data.blob();
        const reader = new FileReader();
        reader.readAsDataURL(blob);
        reader.onloadend = () => {
          this.mutations.setChannelStackImage(reader.result);
        };
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async exportChannelStackImage(format: ExportFormat = "png") {
    if (format === "ome-tiff") {
      try {
        const groupId = this.group?.getters.activeGroupId!;
        const activeAcquisition = this.getters.activeAcquisition;
        const response = await api.downloadOmeTiffImage(groupId, activeAcquisition?.id!);
        const blob = await response.blob();
        const filename = activeAcquisition?.location.replace(/^.*[\\\/]/, "");
        saveAs(blob, filename);
      } catch (error) {
        await this.main!.actions.checkApiError(error);
      }
    } else {
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

  async getAnnotationData() {
    const annotations = this.annotations!.getters.annotations;
    const cellClasses = this.annotations!.getters.cellClasses;
    this.cells!.mutations.updateCellsByAnnotations({ annotations: annotations, cellClasses: cellClasses });
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
      const showMask = this.ui!.getters.showMask;
      const maskOpacity = this.ui!.getters.maskOpacity;
      const activeAcquisitionId = this.getters.activeAcquisitionId;
      if (activeAcquisitionId && activeDataset.meta.masks) {
        const mask = activeDataset.meta.masks[activeAcquisitionId];
        if (mask) {
          output["mask"] = {
            showMask: showMask,
            opacity: maskOpacity,
            location: mask.location,
          };
          if (this.cells?.getters.activeResultId) {
            output["mask"]["resultId"] = this.cells?.getters.activeResultId;
          }

          if (this.cells?.getters.heatmap) {
            output["mask"]["colorsType"] = this.cells.getters.heatmap.type;
            output["mask"]["colorsName"] = this.cells.getters.heatmap.value;
          }

          if (this.cells?.getters.heatmap && this.cells.getters.heatmap.type === "annotation") {
            const annotations = this.annotations?.getters.annotations;
            const cellClasses = this.annotations?.getters.cellClasses;
            const cellColors = {};
            annotations?.forEach((annotation) => {
              if (annotation.visible) {
                const color = cellClasses![annotation.cellClass];
                if (!cellColors.hasOwnProperty(color)) {
                  cellColors[color] = [];
                }
                const objectNumbers: number[] = [];
                annotation.cellIds.forEach((cellId) => {
                  const cell = this.cells?.getters.cells[cellId];
                  if (cell?.acquisitionId === activeAcquisitionId) {
                    objectNumbers.push(cell.objectNumber);
                  }
                });
                cellColors[color] = cellColors[color].concat(objectNumbers);
              }
            });
            output["mask"]["gated"] = true;
            output["mask"]["cells"] = cellColors;
          } else {
            // Prepare selected cell objectNumbers visualisation
            const selectedAcquisitionCells = this.cells?.getters.selectedCells?.filter(
              (v) => v.acquisitionId === activeAcquisitionId
            );
            if (selectedAcquisitionCells && selectedAcquisitionCells.length > 0) {
              output["mask"]["gated"] = true;
              output["mask"]["cells"] = {
                "#ffffff": selectedAcquisitionCells.map((item) => item.objectNumber),
              };
            }
          }
        }
      }
    }
    return output;
  }
}
