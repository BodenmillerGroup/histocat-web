import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import HistogramView from "./HistogramView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class HistogramComponent extends LayoutComponent {
  static readonly typeName = "histogram";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(HistogramView), store, parent);
  }

  protected handleResizeEvent(): void {
    const containerWidth = this._container.width;
    if (containerWidth && containerWidth > 0) {
      this._component.refresh(containerWidth);
    }
  }
}
