import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { UiState } from ".";
import { UiGetters } from "./getters";
import { UiMutations } from "./mutations";

export class UiActions extends Actions<
  UiState,
  UiGetters,
  UiMutations,
  UiActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
