import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import TsneView from "./TsneView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class TsneComponent extends LayoutComponent {
  static readonly typeName = TsneComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(TsneView), store, parent);
  }
}
