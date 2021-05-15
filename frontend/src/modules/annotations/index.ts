import { Module } from "vuex-smart-module";
import { AnnotationsActions } from "./actions";
import { AnnotationsGetters } from "./getters";
import { AnnotationsMutations } from "./mutations";
import { IAnnotation } from "@/modules/annotations/models";

export const defaultCellClasses = {
  none: "#ffffff",
  tumor: "#ff0000",
  stroma: "#37ff00",
  immune: "#b300ff",
  necrosis: "#000000",
  other: "#ffd500",
  region: "#005eff",
  ignore: "#c8c8c8",
  positive: "#ff4800",
  negative: "#c800ff",
};

export class AnnotationsState {
  cellClasses: { [name: string]: string } = defaultCellClasses;
  annotations: IAnnotation[] = [];
}

export const annotationsModule = new Module({
  namespaced: true,

  state: AnnotationsState,
  getters: AnnotationsGetters,
  mutations: AnnotationsMutations,
  actions: AnnotationsActions,
});
