import { Module } from "vuex-smart-module";
import { ResponsiveActions } from "./actions";
import { ResponsiveGetters } from "./getters";
import { ResponsiveMutations } from "./mutations";
import { IResponsive } from "./models";

export class ResponsiveState {
  responsive: IResponsive = {
    width: null,
    height: null,
  };
}

export const responsiveModule = new Module({
  namespaced: true,

  state: ResponsiveState,
  getters: ResponsiveGetters,
  mutations: ResponsiveMutations,
  actions: ResponsiveActions,
});
