import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import RegionView from "./RegionView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class RegionComponent extends LayoutComponent {
  static readonly typeName = "region";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(RegionView), store, parent);
  }
}
