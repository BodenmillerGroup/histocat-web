import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import GatingView from "./GatingView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class GatingComponent extends LayoutComponent {
  static readonly typeName = "gating";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(GatingView), store, parent);
  }
}
