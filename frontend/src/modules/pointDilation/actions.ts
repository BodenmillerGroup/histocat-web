import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { PointDilationState } from ".";
import { PointDilationGetters } from "./getters";
import { PointDilationMutations } from "./mutations";

export class PointDilationActions extends Actions<
  PointDilationState,
  PointDilationGetters,
  PointDilationMutations,
  PointDilationActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
