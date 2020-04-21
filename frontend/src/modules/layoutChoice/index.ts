import { Module } from "vuex-smart-module";
import { LayoutChoiceActions } from "./actions";
import { LayoutChoiceGetters } from "./getters";
import { LayoutChoiceMutations } from "./mutations";

/*
we have a UI heuristic to pick the default layout, based on assumptions
about commonly used names.  Preferentially, pick in the following order:

  1. "umap"
  2. "tsne"
  3. "pca"
  4. give up, use the first available
*/
function bestDefaultLayout(layouts) {
  const preferredNames = ["umap", "tsne", "pca"];
  const idx = preferredNames.findIndex(name => layouts.indexOf(name) !== -1);
  if (idx !== -1) return preferredNames[idx];
  return layouts[0];
}

function setToDefaultLayout(world) {
  const { schema } = world;
  const available = schema.layout.obs.map(v => v.name).sort();
  const current = bestDefaultLayout(available);
  const currentDimNames = schema.layout.obsByName[current].dims;
  return { available, current, currentDimNames };
}

export class LayoutChoiceState {
  available: string[] = []; // all available choices
  current?: string = undefined; // name of the current layout, eg, 'umap'
  currentDimNames: string[] = []; // dimension name
}

export const layoutChoiceModule = new Module({
  namespaced: true,

  state: LayoutChoiceState,
  getters: LayoutChoiceGetters,
  mutations: LayoutChoiceMutations,
  actions: LayoutChoiceActions,
});
