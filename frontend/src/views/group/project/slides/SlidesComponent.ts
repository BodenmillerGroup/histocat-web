import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import SlidesView from "./SlidesView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class SlidesComponent extends LayoutComponent {
  static readonly typeName = SlidesComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(SlidesView), store, parent);
  }
}
