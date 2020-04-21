import { Module } from "vuex-smart-module";
import { ColorsActions } from "./actions";
import { ColorsGetters } from "./getters";
import { ColorsMutations } from "./mutations";

export class ColorsState {
  colorMode: any = null;
  colorAccessor: any = null;
  rgb: any = null;
  scale: any = null;
}

export const colorsModule = new Module({
  namespaced: true,

  state: ColorsState,
  getters: ColorsGetters,
  mutations: ColorsMutations,
  actions: ColorsActions,
});
