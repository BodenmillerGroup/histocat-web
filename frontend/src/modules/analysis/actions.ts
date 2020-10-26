import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnalysisState } from ".";
import { api } from "./api";
import { AnalysisGetters } from "./getters";
import { IRegionStatsSubmission } from "./models";
import { AnalysisMutations } from "./mutations";
import { groupModule } from "@/modules/group";
import { datasetsModule } from "@/modules/datasets";
import { resultsModule } from "@/modules/results";

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  settings?: Context<typeof settingsModule>;
  projects?: Context<typeof projectsModule>;
  datasets?: Context<typeof datasetsModule>;
  results?: Context<typeof resultsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.settings = settingsModule.context(store);
    this.projects = projectsModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.results = resultsModule.context(store);
  }

  async getScatterPlotData(payload: { markerX: string; markerY: string }) {
    try {
      const datasetId = this.datasets!.getters.activeDatasetId!;
      const resultId = this.results!.getters.activeResultId;
      const heatmapType = this.results!.getters.heatmap ? this.results!.getters.heatmap.type : undefined;
      const heatmap = this.results!.getters.heatmap ? this.results!.getters.heatmap.label : undefined;
      const data = await api.getScatterPlotData(
        datasetId,
        resultId,
        payload.markerX,
        payload.markerY,
        heatmapType,
        heatmap
      );
      this.mutations.setScatterPlotData(data);
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

  async getPcaData(resultId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const heatmapType = this.results!.getters.heatmap ? this.results!.getters.heatmap.type : undefined;
      const heatmap = this.results!.getters.heatmap ? this.results!.getters.heatmap.label : undefined;
      const response = await api.getPcaData(groupId, resultId, heatmapType, heatmap);
      this.mutations.setPcaData(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getTsneResult(resultId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const heatmapType = this.results!.getters.heatmap ? this.results!.getters.heatmap.type : undefined;
      const heatmap = this.results!.getters.heatmap ? this.results!.getters.heatmap.label : undefined;
      const data = await api.getTsneData(groupId, resultId, heatmapType, heatmap);
      this.mutations.setTsneData(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getUmapResult(resultId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const heatmapType = this.results!.getters.heatmap ? this.results!.getters.heatmap.type : undefined;
      const heatmap = this.results!.getters.heatmap ? this.results!.getters.heatmap.label : undefined;
      const data = await api.getUmapData(groupId, resultId, heatmapType, heatmap);
      this.mutations.setUmapData(data);
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
