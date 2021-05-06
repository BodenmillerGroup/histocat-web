import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import ScatterView from "./ScatterView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class ScatterComponent extends LayoutComponent {
  static readonly typeName = ScatterComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(ScatterView), store, parent);
  }

  protected handleResizeEvent(): void {
    this._component.refresh();
  }
}
