import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { LayoutChoiceState } from ".";
import { LayoutChoiceGetters } from "./getters";
import { LayoutChoiceMutations } from "./mutations";

export class LayoutChoiceActions extends Actions<
  LayoutChoiceState,
  LayoutChoiceGetters,
  LayoutChoiceMutations,
  LayoutChoiceActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
