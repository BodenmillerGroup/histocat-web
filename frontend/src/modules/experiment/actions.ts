import { datasetModule } from "@/modules/datasets";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { ExperimentState } from ".";
import { api } from "./api";
import { ExperimentGetters } from "./getters";
import { ExportFormat, IExperimentCreate, IExperimentUpdate, IShareCreate } from "./models";
import { ExperimentMutations } from "./mutations";

export class ExperimentActions extends Actions<
  ExperimentState,
  ExperimentGetters,
  ExperimentMutations,
  ExperimentActions
> {
  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;
  datasets?: Context<typeof datasetModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
    this.datasets = datasetModule.context(store);
  }

  async getExperiments() {
    try {
      const data = await api.getExperiments();
      if (data) {
        this.mutations.setExperiments(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getTags() {
    try {
      const data = await api.getTags();
      if (data) {
        this.mutations.setTags(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateExperiment(payload: { id: number; data: IExperimentUpdate }) {
    try {
      const data = await api.updateExperiment(payload.id, payload.data);
      this.mutations.setExperiment(data);
      this.main!.mutations.addNotification({ content: "Experiment successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteExperiment(id: number) {
    try {
      const data = await api.deleteExperiment(id);
      this.mutations.deleteExperiment(id);
      this.main!.mutations.addNotification({ content: "Experiment successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createExperiment(payload: IExperimentCreate) {
    try {
      const data = await api.createExperiment(payload);
      this.mutations.setExperiment(data);
      this.main!.mutations.addNotification({ content: "Experiment successfully created", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async upload(payload: { id: number; data: any }) {
    if (!payload.id) {
      return;
    }
    try {
      await api.upload(
        this.main!.getters.token,
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

  async getExperimentData(id: number) {
    try {
      const data = await api.getExperimentData(id);
      if (data) {
        this.mutations.setExperiment(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStats(payload: { acquisitionId: number; channelName: string }) {
    try {
      return await api.getChannelStats(payload.acquisitionId, payload.channelName);
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
      const response = await api.downloadChannelStackImage(params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        this.mutations.setChannelStackImage(reader.result);
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
      const response = await api.downloadChannelStackImage(params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        this.mutations.setChannelStackImage(reader.result);
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
      const response = await api.downloadChannelStackImage(params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async setSharedChannelLevels(payload: { metal: string; levels: number[] }) {
    const experiment = this.getters.activeExperiment;
    // if (experiment && experiment.slides) {
    //   for (const slide of experiment.slides) {
    //     for (const acquisition of slide.acquisitions) {
    //       for (const channel of acquisition.channels) {
    //         if (channel.name === payload.metal) {
    //           let settings = this.settings!.getters.getChannelSettings(`${channel.acquisition_id}_${channel.name}`);
    //           if (!settings) {
    //             settings = {
    //               id: `${channel.acquisition_id}_${channel.name}`,
    //               customLabel: channel.label,
    //               levels: {
    //                 min: payload.levels[0],
    //                 max: payload.levels[1],
    //               },
    //             };
    //           } else {
    //             settings = {
    //               ...settings,
    //               levels: {
    //                 min: payload.levels[0],
    //                 max: payload.levels[1],
    //               },
    //             };
    //           }
    //           this.settings!.mutations.setChannelSettings(settings);
    //         }
    //       }
    //     }
    //   }
    // }
  }

  async createShare(payload: IShareCreate) {
    try {
      const data = await api.createShare(payload);
      this.main!.mutations.addNotification({ content: "Experiment successfully shared", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getExperimentShares(experimentId: number) {
    try {
      const data = await api.getExperimentShares(experimentId);
      this.mutations.setShares(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteSlide(id: number) {
    try {
      const data = await api.deleteSlide(id);
      this.main!.mutations.addNotification({ content: "Slide successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  private prepareStackParams(format: "png" | "tiff" = "png") {
    const activeAcquisitionId = this.getters.activeAcquisitionId;
    const channels = this.getters.selectedChannels.map((channel) => {
      const color = this.settings!.getters.metalColorMap.get(channel.name);
      const settings = this.settings!.getters.getChannelSettings(activeAcquisitionId, channel.name);
      const min = settings && settings.levels ? settings.levels.min : undefined;
      const max = settings && settings.levels ? settings.levels.max : undefined;
      const customLabel = settings && settings.customLabel ? settings.customLabel : channel.label;
      return {
        name: channel.name,
        color: color,
        customLabel: customLabel,
        min: min,
        max: max,
      };
    });

    const filter = this.settings!.getters.filter;
    const legend = this.settings!.getters.legend;
    const scalebar = this.settings!.getters.scalebar;

    const result = {
      acquisitionId: activeAcquisitionId,
      format: format,
      filter: filter,
      legend: legend,
      scalebar: scalebar,
      channels: channels,
    };

    const activeDataset = this.datasets!.getters.activeDataset;
    if (activeDataset) {
      result["datasetId"] = activeDataset.id;
      const maskSettings = this.settings!.getters.maskSettings;
      const acquisition = this.getters.activeAcquisition;
      if (acquisition && activeDataset.input && activeDataset.input.probability_masks) {
        const mask = activeDataset.input.probability_masks[acquisition.id];
        if (mask) {
          result["mask"] = {
            apply: maskSettings.apply,
            colorize: false,
            location: mask.location,
          };
        }
      }
    }

    return result;
  }
}
