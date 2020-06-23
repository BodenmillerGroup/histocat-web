import { Actions } from "vuex-smart-module";
import { SettingsState } from ".";
import { SettingsGetters } from "./getters";
import { SettingsMutations } from "./mutations";
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_PRESET,
  SET_SHARED_CHANNEL_SETTINGS,
  SET_CHANNEL_SETTINGS,
  SET_METAL_COLOR,
  SET_FILTER,
  SET_LEGEND,
  SET_SCALEBAR,
  SET_MASK_SETTINGS,
} from "./events";
import { IChannelSettings, IImageFilter, IImageLegend, IImageScalebar, IMaskSettings } from "./models";
import { IPreset } from "@/modules/presets/models";

export class SettingsActions extends Actions<SettingsState, SettingsGetters, SettingsMutations, SettingsActions> {
  setSharedChannelSettings(payload: IChannelSettings[], isGlobal = true) {
    BroadcastManager.publish(SET_SHARED_CHANNEL_SETTINGS, payload, isGlobal);
  }

  setChannelSettings(payload: IChannelSettings, isGlobal = true) {
    BroadcastManager.publish(SET_CHANNEL_SETTINGS, payload, isGlobal);
  }

  setMetalColor(payload: { metal: string; color: string }, isGlobal = true) {
    BroadcastManager.publish(SET_METAL_COLOR, payload, isGlobal);
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

  setPreset(payload: IPreset, isGlobal = true) {
    BroadcastManager.publish(SET_PRESET, payload, isGlobal);
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
