import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { ResultsState } from ".";
import { api } from "./api";
import { ResultsGetters } from "./getters";
import { ResultsMutations } from "./mutations";
import { groupModule } from "@/modules/group";

export class ResultsActions extends Actions<ResultsState, ResultsGetters, ResultsMutations, ResultsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
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

  async getResult(resultId: number) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const data = await api.getResult(groupId, resultId);
      this.mutations.setEntity(data);
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
