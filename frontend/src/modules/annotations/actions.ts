import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { AnnotationsState } from ".";
import { AnnotationsGetters } from "./getters";
import { AnnotationsMutations } from "./mutations";
import { datasetsModule } from "@/modules/datasets";
import { groupModule } from "@/modules/group";

export class AnnotationsActions extends Actions<
  AnnotationsState,
  AnnotationsGetters,
  AnnotationsMutations,
  AnnotationsActions
> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;
  dataset?: Context<typeof datasetsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
    this.dataset = datasetsModule.context(store);
  }
}
