import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { GraphSelectionState } from ".";
import { GraphSelectionGetters } from "./getters";
import { GraphSelectionMutations } from "./mutations";

export class GraphSelectionActions extends Actions<GraphSelectionState, GraphSelectionGetters, GraphSelectionMutations, GraphSelectionActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
