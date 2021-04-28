import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnalysisState } from ".";
import { api } from "./api";
import { AnalysisGetters } from "./getters";
import { IClassifyCellsSubmission, IRegionStatsSubmission } from "./models";
import { AnalysisMutations } from "./mutations";
import { groupModule } from "@/modules/group";
import { datasetsModule } from "@/modules/datasets";
import { annotationsModule } from "@/modules/annotations";

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  datasets?: Context<typeof datasetsModule>;
  annotations?: Context<typeof annotationsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.annotations = annotationsModule.context(store);
  }

  async calculateRegionStats(payload: IRegionStatsSubmission) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.calculateRegionStats(groupId, payload);
      this.mutations.setSelectedRegionStats(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async classifyCells(payload: {channels: string[], thresholds: {[cellClass: string]: number}, nEstimators: number}) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const datasetId = this.datasets!.getters.activeDatasetId;
      const cellClasses = this.annotations?.getters.cellClasses;
      const annotations = this.annotations?.getters.annotations;
      if (datasetId && cellClasses && annotations) {
        const params: IClassifyCellsSubmission = {
          dataset_id: datasetId!,
          channels: payload.channels,
          thresholds: payload.thresholds,
          n_estimators: payload.nEstimators,
          cell_classes: cellClasses,
          annotations: annotations,
        };
        const data = await api.classifyCells(groupId, params);
        if (data) {
          this.annotations?.mutations.setCellClasses(data.cellClasses);
          this.annotations?.mutations.setAnnotations(data.annotations);
        }
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
