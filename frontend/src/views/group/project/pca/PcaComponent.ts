import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import PcaView from "./PcaView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class PcaComponent extends LayoutComponent {
  static readonly typeName = "pca";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(PcaView), store, parent);
  }

  protected handleResizeEvent(): void {
    this._component.refresh();
  }
}
