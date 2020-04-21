import { Store } from "vuex";
import { Actions } from "vuex-smart-module";
import { ColorsState } from ".";
import { ColorsGetters } from "./getters";
import { ColorsMutations } from "./mutations";

export class ColorsActions extends Actions<
  ColorsState,
  ColorsGetters,
  ColorsMutations,
  ColorsActions
> {
  // Called after the module is initialized
  $init(store: Store<any>): void {}
}
