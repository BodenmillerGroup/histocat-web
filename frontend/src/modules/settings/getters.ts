import { Getters } from 'vuex-smart-module';
import { SettingsState } from '.';

export class SettingsGetters extends Getters<SettingsState> {
  channelSettings(id: number) {
    return this.state.channelsSettings.get(id);
  }

  get metalColorMap() {
    return this.state.metalColorMap;
  }
}
