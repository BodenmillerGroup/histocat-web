import { Actions } from "vuex-smart-module";
import { SettingsState } from ".";
import { SettingsGetters } from "./getters";
import { SettingsMutations } from "./mutations";
import {BroadcastManager} from "@/utils/BroadcastManager";
import {RESET_SETTINGS} from "@/modules/settings/events";

export class SettingsActions extends Actions<SettingsState, SettingsGetters, SettingsMutations, SettingsActions> {
  async resetSettings() {
    if ("caches" in self) {
      await caches.keys().then((keyList) => {
        return Promise.all(
          keyList.map((key) => {
            return caches.delete(key);
          })
        );
      });
    }
    BroadcastManager.publish(RESET_SETTINGS);
  }
}
