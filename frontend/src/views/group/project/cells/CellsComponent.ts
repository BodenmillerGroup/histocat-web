import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import CellsView from "./CellsView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class CellsComponent extends LayoutComponent {
  static readonly typeName = CellsComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(CellsView), store, parent);
  }
}
