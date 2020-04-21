import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { CentroidLabelsState } from ".";
import { CentroidLabelsGetters } from "./getters";
import { CentroidLabelsMutations } from "./mutations";

export class CentroidLabelsActions extends Actions<
  CentroidLabelsState,
  CentroidLabelsGetters,
  CentroidLabelsMutations,
  CentroidLabelsActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
