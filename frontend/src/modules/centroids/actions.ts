import { Store } from "vuex";
import {Actions, Context} from "vuex-smart-module";
import { CentroidsState } from ".";
import { CentroidsGetters } from "./getters";
import { CentroidsMutations } from "./mutations";
import { api } from "./api";
import { ICentroidsSubmission } from "./models";
import {mainModule} from "@/modules/main";

export class CentroidsActions extends Actions<CentroidsState, CentroidsGetters, CentroidsMutations, CentroidsActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
  }

  async getCentroids(payload: ICentroidsSubmission) {
    try {
      const response = await api.getCentroids(payload);
      this.mutations.setCentroids(response);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
