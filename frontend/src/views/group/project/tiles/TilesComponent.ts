import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import TilesView from "./TilesView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class TilesComponent extends LayoutComponent {
  static readonly typeName = TilesComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(TilesView), store, parent);
  }
}
