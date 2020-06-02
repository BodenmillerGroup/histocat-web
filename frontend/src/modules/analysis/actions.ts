import { experimentModule } from "@/modules/experiment";
import { ExportFormat } from "@/modules/experiment/models";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { saveAs } from "file-saver";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnalysisState } from ".";
import { api } from "./api";
import { AnalysisGetters } from "./getters";
import {
  IImageSegmentationSettings,
  IPCASubmission,
  IPhenoGraphSubmission,
  IRegionStatsSubmission,
  ITSNESubmission,
  IUMAPSubmission,
} from "./models";
import { AnalysisMutations } from "./mutations";

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;
  experiment?: Context<typeof experimentModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
    this.experiment = experimentModule.context(store);
  }

  async getSegmentationImage(settings: IImageSegmentationSettings) {
    const params = await this.actions.prepareSegmentationParams(settings);
    if (params.channels.length === 0) {
      return;
    }
    try {
      const response = await api.produceSegmentationImage(params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        this.mutations.setSegmentationImage(reader.result);
      };
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async exportSegmentationImage(payload: { settings: IImageSegmentationSettings; format: ExportFormat }) {
    const params = await this.actions.prepareSegmentationParams(payload.settings, payload.format);
    try {
      const response = await api.produceSegmentationImage(params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async produceSegmentationContours(settings: IImageSegmentationSettings) {
    const params = await this.actions.prepareSegmentationParams(settings);
    if (params.channels.length === 0) {
      return;
    }
    try {
      const response = await api.produceSegmentationContours(params);
      this.mutations.setSegmentationContours(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getScatterPlotData(payload: {
    datasetId: number;
    acquisitionIds: number[];
    markerX: string;
    markerY: string;
    markerZ: string;
    heatmapType: string;
    heatmap: string;
  }) {
    try {
      const response = await api.getScatterPlotData(
        payload.datasetId,
        payload.acquisitionIds,
        payload.markerX,
        payload.markerY,
        payload.markerZ,
        payload.heatmapType,
        payload.heatmap
      );
      this.mutations.setScatterPlotData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getBoxPlotData(payload: { datasetId: number; acquisitionId: number; markers: string[] }) {
    try {
      const response = await api.getBoxPlotData(payload.datasetId, payload.acquisitionId, payload.markers);
      this.mutations.setBoxPlotData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getPCAData(payload: IPCASubmission) {
    try {
      const response = await api.getPCAData(payload);
      this.mutations.setPCAData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async submitTSNE(payload: ITSNESubmission) {
    try {
      const response = await api.submitTSNE(payload);
      const notification = { content: "t-SNE processing started. This may take a while..." };
      this.main!.mutations.addNotification(notification);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getTSNEResult(payload: { datasetId: number; name: string; heatmapType: string; heatmap: string }) {
    try {
      const response = await api.getTSNEData(payload.datasetId, payload.name, payload.heatmapType, payload.heatmap);
      this.mutations.setTSNEData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async submitUMAP(payload: IUMAPSubmission) {
    try {
      const response = await api.submitUMAP(payload);
      const notification = { content: "UMAP processing started. This may take a while..." };
      this.main!.mutations.addNotification(notification);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getUMAPResult(payload: { datasetId: number; name: string; heatmapType: string; heatmap: string }) {
    try {
      const response = await api.getUMAPData(payload.datasetId, payload.name, payload.heatmapType, payload.heatmap);
      this.mutations.setUMAPData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async submitPhenoGraph(payload: IPhenoGraphSubmission) {
    try {
      const response = await api.submitPhenoGraph(payload);
      const notification = { content: "PhenoGraph processing started. This may take a while..." };
      this.main!.mutations.addNotification(notification);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getPhenoGraphResult(payload: { datasetId: number; name: string }) {
    try {
      const response = await api.getPhenoGraphData(payload.datasetId, payload.name);
      this.mutations.setPhenoGraphData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async calculateRegionStats(payload: IRegionStatsSubmission) {
    try {
      const response = await api.calculateRegionStats(payload);
      this.mutations.setSelectedRegionStats(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  private prepareSegmentationParams(segmentationSettings: IImageSegmentationSettings, format: "png" | "tiff" = "png") {
    const activeAcquisitionId = this.experiment!.getters.activeAcquisitionId;
    const channels = this.experiment!.getters.selectedChannels.map((channel) => {
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
    const scalebar = this.settings!.getters.scalebar;

    return {
      format: format,
      filter: filter,
      scalebar: scalebar,
      channels: channels,
      settings: segmentationSettings,
    };
  }
}
