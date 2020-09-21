import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnalysisState } from ".";
import { api } from "./api";
import { AnalysisGetters } from "./getters";
import {
  IPCASubmission,
  IPhenoGraphSubmission,
  IRegionStatsSubmission,
  ITSNESubmission,
  IUMAPSubmission,
} from "./models";
import { AnalysisMutations } from "./mutations";
import { groupModule } from "@/modules/group";

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  projects?: Context<typeof projectsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.projects = projectsModule.context(store);
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

  async getBoxPlotData(payload: {
    datasetId: number;
    gateId: number | null;
    acquisitionIds: number[];
    markers: string[];
  }) {
    try {
      const response = await api.getBoxPlotData(
        payload.datasetId,
        payload.gateId,
        payload.acquisitionIds,
        payload.markers
      );
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

  async getTSNEResult(payload: { resultId: number; heatmapType: string; heatmap: string }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getTSNEData(groupId, payload.resultId, payload.heatmapType, payload.heatmap);
      this.mutations.setTSNEData(data);
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

  async getUMAPResult(payload: { resultId: number; heatmapType: string; heatmap: string }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getUMAPData(groupId, payload.resultId, payload.heatmapType, payload.heatmap);
      this.mutations.setUMAPData(data);
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

  async getPhenoGraphResult(payload: { resultId: number }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getPhenoGraphData(groupId, payload.resultId);
      this.mutations.setPhenoGraphData(data);
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
}
