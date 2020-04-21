import { Module } from "vuex-smart-module";
import { ControlsActions } from "./actions";
import { ControlsGetters } from "./getters";
import { ControlsMutations } from "./mutations";

export class ControlsState {
  // data loading flag
  loading = true;
  error: any = null;

  // all of the data + selection state
  userDefinedGenes: any[] = [];
  userDefinedGenesLoading = false;
  diffexpGenes: any[] = [];

  resettingInterface = false;
  graphInteractionMode = "select";
  opacityForDeselectedCells = 0.2;

  // just easier to read
  scatterplotXXaccessor: any = null;
  scatterplotYYaccessor: any = null;

  // integer as <Component key={graphRenderCounter} - a change in key forces a remount
  graphRenderCounter = 0;
}

export const controlsModule = new Module({
  namespaced: true,

  state: ControlsState,
  getters: ControlsGetters,
  mutations: ControlsMutations,
  actions: ControlsActions,
});
