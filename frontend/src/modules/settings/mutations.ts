import { Mutations } from 'vuex-smart-module';
import { SettingsState } from '.';
import { IChannelSettings } from './models';


export class SettingsMutations extends Mutations<SettingsState> {
  setChannelSettings(channelSettings: IChannelSettings) {
    this.state.channelsSettings.set(channelSettings.id, channelSettings);
  }

  setMetalColor(payload: { metal: string, color: string }) {
    this.state.metalColorMap.set(payload.metal, payload.color);
    this.state.metalColorMap = new Map(this.state.metalColorMap);
  }
}
