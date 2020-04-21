import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { UniverseState } from ".";
import { UniverseGetters } from "./getters";
import { UniverseMutations } from "./mutations";

export class UniverseActions extends Actions<UniverseState, UniverseGetters, UniverseMutations, UniverseActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
