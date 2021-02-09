import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnalysisState } from ".";
import { api } from "./api";
import { AnalysisGetters } from "./getters";
import { IRegionStatsSubmission } from "./models";
import { AnalysisMutations } from "./mutations";
import { groupModule } from "@/modules/group";

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
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
}
