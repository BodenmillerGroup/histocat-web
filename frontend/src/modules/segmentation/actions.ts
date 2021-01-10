import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { SegmentationState } from ".";
import { api } from "./api";
import { SegmentationGetters } from "./getters";
import { SegmentationMutations } from "./mutations";
import { groupModule } from "@/modules/group";

export class SegmentationActions extends Actions<
  SegmentationState,
  SegmentationGetters,
  SegmentationMutations,
  SegmentationActions
> {
  // Declare context type
  main?: Context<typeof mainModule>;
  group?: Context<typeof groupModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.group = groupModule.context(store);
  }
}
