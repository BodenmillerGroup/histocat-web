import { Module } from "vuex-smart-module";
import { PresetActions } from "./actions";
import { PresetGetters } from "./getters";
import { IPreset } from "./models";
import { PresetMutations } from "./mutations";
import { schema } from "normalizr";

export const presetSchema = new schema.Entity("presets");
export const presetListSchema = [presetSchema];

export class PresetState {
  ids: ReadonlyArray<number> = [];
  entities: { [key: number]: IPreset } = {};
  activePresetId: number | null = null;
}

export const presetModule = new Module({
  namespaced: true,

  state: PresetState,
  getters: PresetGetters,
  mutations: PresetMutations,
  actions: PresetActions,
});
