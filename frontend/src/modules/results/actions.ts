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
import { IResultUpdate, ISelectedCell } from "@/modules/results/models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SELECTED_CELLS } from "@/modules/results/events";

export class ResultsActions extends Actions<ResultsState, ResultsGetters, ResultsMutations, ResultsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  pipelines?: Context<typeof pipelinesModule>;
  analysis?: Context<typeof analysisModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.pipelines = pipelinesModule.context(store);
    this.analysis = analysisModule.context(store);
  }

  setSelectedCells(payload: ISelectedCell[], isGlobal = true) {
    BroadcastManager.publish(SET_SELECTED_CELLS, payload, isGlobal);
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

  async loadResultData(resultId: number) {
    try {
      this.mutations.setActiveResultId(resultId);
      this.analysis?.mutations.reset();
      const groupId = this.group?.getters.activeGroupId!;
      const colorsType = this.getters.heatmap ? this.getters.heatmap.type : undefined;
      const colorsName = this.getters.heatmap ? this.getters.heatmap.label : undefined;
      const data = await api.getResultData(groupId, resultId, colorsType, colorsName);
      this.mutations.setCells(data);
      const result = this.getters.activeResult;
      if (result) {
        this.pipelines?.mutations.setSteps(result.pipeline);
        this.pipelines?.mutations.setSelectedAcquisitionIds(result.input);
      }
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
