import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import PresetsView from "./PresetsView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class PresetsComponent extends LayoutComponent {
  static readonly typeName = "presets";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(PresetsView), store, parent);
  }
}
