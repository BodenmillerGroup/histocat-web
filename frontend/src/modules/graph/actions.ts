import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GraphState } from ".";
import { GraphGetters } from "./getters";
import { GraphMutations } from "./mutations";

export class GraphActions extends Actions<GraphState, GraphGetters, GraphMutations, GraphActions> {
  // Declare context type
  main?: Context<typeof mainModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
  }
}
