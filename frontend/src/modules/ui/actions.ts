import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { UiState } from ".";
import { UiGetters } from "./getters";
import { UiMutations } from "./mutations";
import { settingsModule } from "@/modules/settings";

export class UiActions extends Actions<UiState, UiGetters, UiMutations, UiActions> {
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.settings = settingsModule.context(store);
  }

  addLayout(name: string) {
    this.mutations.addLayout(name);
    this.settings?.mutations.setActiveLayoutUid(this.getters.activeLayout.uid);
  }

  resetLayouts() {
    this.mutations.resetLayouts();
    this.settings?.mutations.setActiveLayoutUid(this.getters.activeLayout.uid);
  }

  loadLayout(uid: string) {
    this.mutations.loadLayout(uid);
    this.settings?.mutations.setActiveLayoutUid(this.getters.activeLayout.uid);
  }
}
