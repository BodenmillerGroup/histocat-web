import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { ControlsState } from ".";
import { ControlsGetters } from "./getters";
import { ControlsMutations } from "./mutations";

export class ControlsActions extends Actions<ControlsState, ControlsGetters, ControlsMutations, ControlsActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
