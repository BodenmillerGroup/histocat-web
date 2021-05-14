import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import LeidenView from "./LeidenView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class LeidenComponent extends LayoutComponent {
  static readonly typeName = "leiden";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(LeidenView), store, parent);
  }
}
