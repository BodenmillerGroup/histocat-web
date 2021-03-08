import { Actions } from "vuex-smart-module";
import { SettingsState } from ".";
import { SettingsGetters } from "./getters";
import { SettingsMutations } from "./mutations";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_CHANNELS_SETTINGS, SET_FILTER, SET_LEGEND, SET_SCALEBAR, SET_MASK_SETTINGS } from "./events";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";

export class SettingsActions extends Actions<SettingsState, SettingsGetters, SettingsMutations, SettingsActions> {
  setChannelsSettings(payload: { channelName: string; settings: IChannelSettings }, isGlobal = true) {
    BroadcastManager.publish(SET_CHANNELS_SETTINGS, payload, isGlobal);
  }

  setFilter(payload: IImageFilter, isGlobal = true) {
    BroadcastManager.publish(SET_FILTER, payload, isGlobal);
  }

  setLegend(payload: IImageLegend, isGlobal = true) {
    BroadcastManager.publish(SET_LEGEND, payload, isGlobal);
  }

  setScalebar(payload: IImageScalebar, isGlobal = true) {
    BroadcastManager.publish(SET_SCALEBAR, payload, isGlobal);
  }

  setMaskSettings(payload: IMaskSettings, isGlobal = true) {
    BroadcastManager.publish(SET_MASK_SETTINGS, payload, isGlobal);
  }

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
    this.mutations.reset();
  }
}
