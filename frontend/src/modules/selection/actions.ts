import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { SelectionState } from ".";
import { SelectionGetters } from "./getters";
import { SelectionMutations } from "./mutations";

export class SelectionActions extends Actions<
  SelectionState,
  SelectionGetters,
  SelectionMutations,
  SelectionActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
