import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import SettingsView from "./SettingsView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class SettingsComponent extends LayoutComponent {
  static readonly typeName = SettingsComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(SettingsView), store, parent);
  }
}
