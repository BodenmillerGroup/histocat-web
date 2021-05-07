import { Getters } from "vuex-smart-module";
import { SettingsState } from ".";

export class SettingsGetters extends Getters<SettingsState> {
  get activeLayoutUid() {
    return this.state.activeLayoutUid;
  }

  get channelsSettings() {
    return this.state.channelsSettings;
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
}
