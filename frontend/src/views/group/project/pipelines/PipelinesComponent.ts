import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import PipelinesView from "./PipelinesView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class PipelinesComponent extends LayoutComponent {
  static readonly typeName = PipelinesComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(PipelinesView), store, parent);
  }
}
