import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { ResultsState } from ".";
import { api } from "./api";
import { ResultsGetters } from "./getters";
import { ResultsMutations } from "./mutations";
import { groupModule } from "@/modules/group";
import { pipelinesModule } from "@/modules/pipelines";
import { analysisModule } from "@/modules/analysis";
import { IResultUpdate } from "@/modules/results/models";
import { datasetsModule } from "@/modules/datasets";

export class ResultsActions extends Actions<ResultsState, ResultsGetters, ResultsMutations, ResultsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  pipelines?: Context<typeof pipelinesModule>;
  analysis?: Context<typeof analysisModule>;
  datasets?: Context<typeof datasetsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.pipelines = pipelinesModule.context(store);
    this.analysis = analysisModule.context(store);
    this.datasets = datasetsModule.context(store);
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
      this.mutations.setActiveResultId(resultId);
      this.mutations.resetResultData();
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getResultData(groupId, resultId);
      this.mutations.setData(data);
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
        this.mutations.setColors(null);
        return;
      }
      const groupId = this.group?.getters.activeGroupId!;
      const resultId = this.getters.activeResultId!;
      const data = await api.getColorsData(groupId, resultId, colorsType, colorsName);
      this.mutations.setColors(data);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getScatterPlotData(payload: { markerX: string; markerY: string }) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const resultId = this.getters.activeResultId!;
      const data = await api.getScatterPlotData(groupId, resultId, payload.markerX, payload.markerY);
      this.mutations.setScatterData(data);
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
