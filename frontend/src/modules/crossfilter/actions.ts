import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { CrossfilterState } from ".";
import { CrossfilterGetters } from "./getters";
import { CrossfilterMutations } from "./mutations";

export class CrossfilterActions extends Actions<CrossfilterState, CrossfilterGetters, CrossfilterMutations, CrossfilterActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
