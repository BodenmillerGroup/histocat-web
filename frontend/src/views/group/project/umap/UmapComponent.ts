import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import UmapView from "./UmapView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class UmapComponent extends LayoutComponent {
  static readonly typeName = "umap";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(UmapView), store, parent);
  }

  protected handleResizeEvent(): void {
    this._component.refresh();
  }
}
