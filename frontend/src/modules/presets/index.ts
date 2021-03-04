import { Module } from "vuex-smart-module";
import { PresetsActions } from "./actions";
import { PresetsGetters } from "./getters";
import { IPreset } from "./models";
import { PresetsMutations } from "./mutations";
import { schema } from "normalizr";

export const presetSchema = new schema.Entity("presets");
export const presetListSchema = [presetSchema];

export class PresetsState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IPreset } = {};
  activePresetId: number | null = null;
}

export const presetsModule = new Module({
  namespaced: true,

  state: PresetsState,
  getters: PresetsGetters,
  mutations: PresetsMutations,
  actions: PresetsActions,
});
