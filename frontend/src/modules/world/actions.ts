import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { WorldState } from ".";
import { WorldGetters } from "./getters";
import { WorldMutations } from "./mutations";

export class WorldActions extends Actions<WorldState, WorldGetters, WorldMutations, WorldActions> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
