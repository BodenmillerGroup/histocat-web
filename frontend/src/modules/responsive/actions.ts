import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { ResponsiveState } from ".";
import { ResponsiveGetters } from "./getters";
import { ResponsiveMutations } from "./mutations";

export class ResponsiveActions extends Actions<
  ResponsiveState,
  ResponsiveGetters,
  ResponsiveMutations,
  ResponsiveActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
