import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { CentroidsState } from ".";
import { CentroidsGetters } from "./getters";
import { CentroidsMutations } from "./mutations";
import { api } from "./api";
import { ICentroidsData, ICentroidsSubmission } from "./models";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_CENTROIDS } from "./events";
import { groupModule } from "@/modules/group";

export class CentroidsActions extends Actions<CentroidsState, CentroidsGetters, CentroidsMutations, CentroidsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
  }

  setCentroids(payload: ICentroidsData, isGlobal = true) {
    BroadcastManager.publish(SET_CENTROIDS, payload, isGlobal);
  }

  async getCentroids(payload: ICentroidsSubmission) {
    try {
      const groupId = this.group?.getters.activeGroupId!;
      const response = await api.getCentroids(groupId, payload);
      this.actions.setCentroids(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
