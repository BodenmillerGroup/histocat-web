import create from "zustand";
import { schema, normalize } from "normalizr";
import {
  ExportFormat,
  IAcquisition,
  IChannel,
  IChannelStats,
  IChannelUpdate,
  IProject,
  IProjectCreate,
  IProjectData,
  IProjectUpdate,
} from "./models";
import { api } from "./api";
import { displayApiError } from "utils/api";
import { isEqual } from "lodash-es";
import { useGroupsStore } from "../groups";
import { AppToaster } from "../../utils/toaster";
import { useAuthStore } from "../auth";
import { useSettingsStore } from "../settings";
import { useDatasetsStore } from "../datasets";
import { useResultsStore } from "../results";
import { Classes, IToastProps, ProgressBar } from "@blueprintjs/core";

const projectSchema = new schema.Entity("projects");
const projectListSchema = [projectSchema];

type ProjectsState = {
  ids: ReadonlyArray<number>;
  entities: { [key: number]: IProject };
  projectsTags: string[];
  projectData: IProjectData | null;
  activeProjectId: number | null;
  activeAcquisitionId: number | null;
  activeWorkspaceNode: object | null;
  selectedMetals: string[];
  channelStackImage: string | ArrayBuffer | null;
  colorizeMaskInProgress: boolean;

  setActiveProjectId(id: number | null): void;
  setActiveAcquisitionId(id: number | null): void;
  setActiveWorkspaceNode(node: { id: number; type: string } | null): void;
  setSelectedMetals(metals: string[]): void;
  getGroupProjects(groupId: number): Promise<void>;
  getProjectsTags(groupId: number): Promise<void>;
  updateProject(id: number, params: IProjectUpdate): Promise<void>;
  deleteProject(projectId: number): Promise<void>;
  createProject(params: IProjectCreate): Promise<void>;
  uploadSlide(data: FormData): Promise<void>;
  getProject(projectId: number): Promise<void>;
  getProjectData(projectId: number): Promise<void>;
  getChannelStats(acquisitionId: number, channelName: string): Promise<IChannelStats | undefined>;
  getChannelStackImage(): Promise<void>;
  exportChannelStackImage(format: ExportFormat): Promise<void>;
  deleteSlide(id: number): Promise<void>;
  updateChannel(acquisitionId: number, params: IChannelUpdate): Promise<void>;
  getActiveAcquisition(): IAcquisition | undefined;
  getSelectedChannels(): IChannel[];
  prepareStackParams(format?: "png" | "tiff"): any;
};

export const useProjectsStore = create<ProjectsState>((set, get) => ({
  ids: [],
  entities: {},
  activeProjectId: null,
  projectsTags: [],
  projectData: null,
  activeAcquisitionId: null,
  activeWorkspaceNode: null,
  selectedMetals: [],
  channelStackImage: null,
  colorizeMaskInProgress: false,

  setActiveProjectId(id: number | null) {
    set({ activeProjectId: id });
  },

  setActiveAcquisitionId(id: number | null) {
    set({ activeAcquisitionId: id });
  },

  setActiveWorkspaceNode(node: { id: number; type: string } | null) {
    set({ activeWorkspaceNode: node });
  },

  setSelectedMetals(metals: string[]) {
    set({ selectedMetals: metals });
  },

  async getGroupProjects(groupId: number) {
    try {
      const data = await api.getGroupProjects(groupId);
      if (data) {
        const normalizedData = normalize<IProject>(data, projectListSchema);
        set({
          ids: normalizedData.result,
          entities: normalizedData.entities.projects ? Object.freeze(normalizedData.entities.projects) : {},
        });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getProjectsTags(groupId: number) {
    try {
      const data = await api.getProjectsTags(groupId);
      if (data) {
        if (!isEqual(data, get().projectsTags)) {
          set({ projectsTags: data });
        }
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateProject(id: number, params: IProjectUpdate) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.updateProject(groupId, id, params);
      if (data) {
        set({ entities: Object.freeze({ ...get().entities, [data.id]: data }) });
        AppToaster.show({ message: "Project successfully updated", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteProject(projectId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.deleteProject(groupId, projectId);
      if (data) {
        const entities = { ...get().entities };
        delete entities[projectId];
        set({ ids: get().ids.filter((item) => item !== projectId), entities: Object.freeze(entities) });
        AppToaster.show({ message: "Project successfully deleted", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async createProject(params: IProjectCreate) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.createProject(groupId, params);
      if (data) {
        set({ ids: get().ids.concat(data.id), entities: Object.freeze({ ...get().entities, [data.id]: data }) });
        AppToaster.show({ message: "Project successfully created", intent: "success" });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async uploadSlide(data: FormData) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const activeProjectId = get().activeProjectId!;
      const token = useAuthStore.getState().token!;
      AppToaster.show({ message: `File upload: ${0}`, intent: "success" }, "slideUploadToast");
      const xhr = api.uploadSlide(
        token,
        groupId,
        activeProjectId,
        data,
        () => {
          const toast: IToastProps = {
            className: Classes.DARK,
            icon: "cloud-upload",
            message: <ProgressBar value={0} />,
            onDismiss: (didTimeoutExpire: boolean) => {
              if (!didTimeoutExpire) {
                // user dismissed toast with click
                xhr.abort();
              }
            },
            timeout: 0,
          };
          AppToaster.show(toast, "slideUploadToast");
        },
        () => {
          const toast: IToastProps = {
            icon: "cloud-upload",
            message: "Slide successfully uploaded",
            intent: "success",
            timeout: 2000,
          };
          AppToaster.show(toast, "slideUploadToast");
        },
        (event) => {
          const percent = Math.round((100 * event.loaded) / event.total);
          const toast: IToastProps = {
            className: Classes.DARK,
            icon: "cloud-upload",
            message: <ProgressBar value={percent / 100} />,
            onDismiss: (didTimeoutExpire: boolean) => {
              if (!didTimeoutExpire) {
                // user dismissed toast with click
                xhr.abort();
              }
            },
            timeout: percent < 100 ? 0 : 2000,
          };
          AppToaster.show(toast, "slideUploadToast");
        },
        () => {}
      );
    } catch (error) {
      displayApiError(error);
    }
  },

  async getProject(projectId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getProject(groupId, projectId);
      if (data) {
        const existingId = get().ids.find((id) => id === data.id);
        if (!existingId) {
          set({ ids: get().ids.concat(data.id) });
        }
        set({ entities: Object.freeze({ ...get().entities, [data.id]: data }) });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getProjectData(projectId: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.getProjectData(groupId, projectId);
      if (data) {
        set({ projectData: data });
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async getChannelStats(acquisitionId: number, channelName: string) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      return await api.getChannelStats(groupId, acquisitionId, channelName);
    } catch (error) {
      displayApiError(error);
    }
  },

  async getChannelStackImage() {
    const params = get().prepareStackParams();
    if (params.channels.length === 0) {
      return;
    }
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const response = await api.downloadChannelStackImage(groupId, params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        set({ channelStackImage: reader.result });
      };
    } catch (error) {
      displayApiError(error);
    }
  },

  async exportChannelStackImage(format: ExportFormat = "png") {
    const params = get().prepareStackParams(format);
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const response = await api.downloadChannelStackImage(groupId, params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      displayApiError(error);
    }
  },

  async deleteSlide(id: number) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.deleteSlide(groupId, id);
      if (data) {
        AppToaster.show({ message: "Slide successfully deleted", intent: "success" });
        const activeProjectId = get().activeProjectId!;
        await get().getProjectData(activeProjectId);
      }
    } catch (error) {
      displayApiError(error);
    }
  },

  async updateChannel(acquisitionId: number, params: IChannelUpdate) {
    try {
      const groupId = useGroupsStore.getState().activeGroupId!;
      const data = await api.updateChannel(groupId, acquisitionId, params);
      if (data) {
        set({ projectData: data });
      }
      await get().getChannelStackImage();
      AppToaster.show({ message: "Channel successfully updated", intent: "success" });
    } catch (error) {
      displayApiError(error);
    }
  },

  getActiveAcquisition() {
    const projectData = get().projectData;
    const activeAcquisitionId = get().activeAcquisitionId;
    if (projectData && projectData.slides) {
      for (const slide of projectData.slides) {
        const acquisition = slide.acquisitions.find((item) => item.id === activeAcquisitionId);
        if (acquisition) {
          return acquisition;
        }
      }
    }
    return undefined;
  },

  getSelectedChannels() {
    const activeAcquisition = get().getActiveAcquisition();
    const selectedMetals = get().selectedMetals;
    if (activeAcquisition) {
      return Object.values(activeAcquisition.channels).filter((channel) => selectedMetals.includes(channel.name));
    } else {
      return [];
    }
  },

  /**
   * Prepare stack image call parameters
   */
  prepareStackParams(format: "png" | "tiff" = "png") {
    const activeAcquisitionId = get().activeAcquisitionId;
    const { channelsSettings, filter, scalebar } = useSettingsStore.getState();
    const channels = get()
      .getSelectedChannels()
      .map((channel) => {
        const channelSettings = channelsSettings[channel.name];
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

    const output: any = {
      acquisitionId: activeAcquisitionId,
      format: format,
      filter: filter,
      scalebar: scalebar,
      channels: channels,
    };

    const activeDataset = useDatasetsStore.getState().getActiveDataset();
    if (activeDataset) {
      output["datasetId"] = activeDataset.id;
      const maskSettings = useSettingsStore.getState().mask;
      if (activeAcquisitionId && activeDataset.meta.masks) {
        const mask = activeDataset.meta.masks[activeAcquisitionId];
        if (mask) {
          output["mask"] = {
            colorize: false,
            mode: maskSettings.mode,
            location: mask.location,
          };
          const { heatmap, activeResultId, selectedCells } = useResultsStore.getState();
          if (heatmap) {
            output["mask"]["colorsType"] = heatmap.type;
            output["mask"]["colorsName"] = heatmap.label;
          }
          if (activeResultId) {
            output["mask"]["resultId"] = activeResultId;
          }
          // Prepare selected cell ids visualisation
          const selectedCellsFiltered = selectedCells?.filter((v) => v.acquisitionId === activeAcquisitionId);
          if (selectedCellsFiltered && selectedCellsFiltered.length > 0) {
            output["mask"]["gated"] = true;
            output["mask"]["cellIds"] = selectedCellsFiltered.map((item) => item.objectNumber);
          }
        }
      }
    }
    return output;
  },
}));
