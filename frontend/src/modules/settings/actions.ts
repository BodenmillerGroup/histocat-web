import { Actions } from 'vuex-smart-module';
import { SettingsState } from '.';
import { SettingsGetters } from './getters';
import { SettingsMutations } from './mutations';

export class SettingsActions extends Actions<SettingsState, SettingsGetters, SettingsMutations, SettingsActions> {
}
