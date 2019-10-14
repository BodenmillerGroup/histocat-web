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
      const data = await api.getExperiments(this.main!.getters.token);
      if (data) {
        this.mutations.setExperiments(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getTags() {
    try {
      const data = await api.getTags(this.main!.getters.token);
      if (data) {
        this.mutations.setTags(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateExperiment(payload: { id: number; data: IExperimentUpdate }) {
    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.updateExperiment(this.main!.getters.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setExperiment(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Experiment successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteExperiment(id: number) {
    try {
      const notification = { content: "deleting", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteExperiment(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.deleteExperiment(id);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Experiment successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async createExperiment(payload: IExperimentCreate) {
    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createExperiment(this.main!.getters.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.mutations.setExperiment(data);
      this.main!.mutations.removeNotification(notification);
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
        event => {
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
      const data = await api.getExperimentData(this.main!.getters.token, id);
      if (data) {
        this.mutations.setExperiment(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStats(id: number) {
    try {
      return await api.getChannelStats(this.main!.getters.token, id);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getChannelStackImage() {
    const params = this.prepareStackParams();
    if (params.channels.length === 0) {
      return;
    }
    try {
      const response = await api.downloadChannelStackImage(this.main!.getters.token, params);
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
    const params = this.prepareStackParams();
    if (params.channels.length === 0 || !params.hasOwnProperty("mask")) {
      return;
    }
    params["mask"]["apply"] = true;
    params["mask"]["colorize"] = true;
    try {
      this.mutations.setColorizeMaskInProgress(true);
      const response = await api.downloadChannelStackImage(this.main!.getters.token, params);
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
    const params = this.prepareStackParams(format);
    try {
      const response = await api.downloadChannelStackImage(this.main!.getters.token, params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async setSharedChannelLevels(payload: { metal: string; levels: number[] }) {
    const experiment = this.getters.activeExperiment;
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            for (const acquisition of roi.acquisitions) {
              for (const channel of acquisition.channels) {
                if (channel.metal === payload.metal) {
                  let settings = this.settings!.getters.getChannelSettings(channel.id);
                  if (!settings) {
                    settings = {
                      id: channel.id,
                      customLabel: channel.label,
                      levels: {
                        min: payload.levels[0],
                        max: payload.levels[1]
                      }
                    };
                  } else {
                    settings = {
                      ...settings,
                      levels: {
                        min: payload.levels[0],
                        max: payload.levels[1]
                      }
                    };
                  }
                  this.settings!.mutations.setChannelSettings(settings);
                }
              }
            }
          }
        }
      }
    }
  }

  async createShare(payload: IShareCreate) {
    try {
      const notification = { content: "saving", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createShare(this.main!.getters.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Experiment successfully shared", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getExperimentShares(experimentId: number) {
    try {
      const data = await api.getExperimentShares(this.main!.getters.token, experimentId);
      this.mutations.setShares(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteSlide(id: number) {
    try {
      const notification = { content: "deleting", showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteSlide(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500))
      ]))[0];
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: "Slide successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  private prepareStackParams(format: "png" | "tiff" = "png") {
    const channels = this.getters.selectedChannels.map(channel => {
      const color = this.settings!.getters.metalColorMap.get(channel.metal);
      const settings = this.settings!.getters.getChannelSettings(channel.id);
      const min = settings && settings.levels ? settings.levels.min : undefined;
      const max = settings && settings.levels ? settings.levels.max : undefined;
      const customLabel = settings && settings.customLabel ? settings.customLabel : channel.label;
      return {
        id: channel.id,
        color: color,
        customLabel: customLabel,
        min: min,
        max: max
      };
    });

    const filter = this.settings!.getters.filter;
    const legend = this.settings!.getters.legend;
    const scalebar = this.settings!.getters.scalebar;

    const result = {
      format: format,
      filter: filter,
      legend: legend,
      scalebar: scalebar,
      channels: channels
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
            location: mask.location
          };
        }
      }
    }

    return result;
  }
}
