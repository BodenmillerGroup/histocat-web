import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { CellsState } from ".";
import { CellsGetters } from "./getters";
import { CellsMutations } from "./mutations";
import { api } from "./api";
import { ICentroidsSubmission, IResultUpdate } from "./models";
import { mainModule } from "@/modules/main";
import { groupModule } from "@/modules/group";
import { pipelinesModule } from "@/modules/pipelines";
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { annotationsModule } from "@/modules/annotations";

export class CellsActions extends Actions<CellsState, CellsGetters, CellsMutations, CellsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  pipelines?: Context<typeof pipelinesModule>;
  analysis?: Context<typeof analysisModule>;
  datasets?: Context<typeof datasetsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.pipelines = pipelinesModule.context(store);
    this.analysis = analysisModule.context(store);
    this.datasets = datasetsModule.context(store);
  }

  async initializeCells(payload: ICentroidsSubmission) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.getCentroids(groupId, payload);
      this.mutations.initializeCells(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getDatasetResults(datasetId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getDatasetResults(groupId, datasetId);
      this.mutations.setEntities(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getResultData(resultId: number) {
    try {
      this.mutations.resetResultData();
      this.mutations.setActiveResultId(resultId);
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getResultData(groupId, resultId);
      this.mutations.updateCellsByResult(data);
      const result = this.getters.activeResult;
      if (result) {
        this.pipelines?.mutations.setSteps(result.pipeline);
        this.pipelines?.mutations.setSelectedAcquisitionIds(result.input);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getColorsData() {
    try {
      const colorsType = this.getters.heatmap ? this.getters.heatmap.type : undefined;
      const colorsName = this.getters.heatmap ? this.getters.heatmap.value : undefined;
      if (!colorsType || !colorsName) {
        this.mutations.updateCellsByColors(null);
        return;
      }
      const groupId = this.group?.getters.activeGroupId!;
      const resultId = this.getters.activeResultId;
      if (resultId) {
        const data = await api.getColorsData(groupId, resultId, colorsType, colorsName);
        this.mutations.updateCellsByColors(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getScatterPlotData(payload: { markerX: string; markerY: string }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const resultId = this.getters.activeResultId!;
      const data = await api.getScatterPlotData(groupId, resultId, payload.markerX, payload.markerY);
      this.mutations.updateCellsByScatterplot(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateResult(payload: { resultId: number; data: IResultUpdate }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.updateResult(groupId, payload.resultId, payload.data);
      this.mutations.updateEntity(data);
      this.main!.mutations.addNotification({ content: "Result successfully updated", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteResult(resultId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.deleteResult(groupId, resultId);
      this.mutations.deleteEntity(resultId);
      this.main!.mutations.addNotification({ content: "Result successfully deleted", color: "success" });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
