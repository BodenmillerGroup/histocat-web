import { Getters } from "vuex-smart-module";
import { SettingsState } from ".";

export class SettingsGetters extends Getters<SettingsState> {
  getChannelSettings(acquisitionId: number, channelName: string) {
    const acquisitionSettings = this.state.channelSettings[acquisitionId];
    return acquisitionSettings ? acquisitionSettings[channelName] : undefined;
  }

  get channelSettings() {
    return this.state.channelSettings;
  }

  get colorMap() {
    return this.state.colorMap;
  }

  get filter() {
    return this.state.filter;
  }

  get legend() {
    return this.state.legend;
  }

  get scalebar() {
    return this.state.scalebar;
  }

  get segmentationSettings() {
    return this.state.segmentationSettings;
  }

  get maskSettings() {
    return this.state.mask;
  }
}
