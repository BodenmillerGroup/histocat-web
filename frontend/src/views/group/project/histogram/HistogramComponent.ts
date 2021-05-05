import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import HistogramView from "./HistogramView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class HistogramComponent extends LayoutComponent {
  static readonly typeName = HistogramComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(HistogramView), store, parent);
  }

  protected handleResizeEvent(): void {
    const containerWidth = this._container.width;
    this._component.refresh(containerWidth);
  }
}
