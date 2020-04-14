import { Getters } from "vuex-smart-module";
import { SettingsState } from ".";

export class SettingsGetters extends Getters<SettingsState> {
  getChannelSettings(acquisitionId: number | undefined, channelName: string) {
    if (!acquisitionId) {
      return undefined;
    }
    const acquisitionSettings = this.state.channelSettings.get(acquisitionId);
    return acquisitionSettings ? acquisitionSettings.get(channelName) : undefined;
  }

  get metalColorMap() {
    return this.state.metalColorMap;
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
