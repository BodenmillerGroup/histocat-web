import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import ImageView from "./ImageView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class ImageComponent extends LayoutComponent {
  static readonly typeName = "image";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(ImageView), store, parent);
  }

  protected handleResizeEvent(): void {
    this._component.refresh();
  }
}
