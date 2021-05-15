import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import SegmentationView from "./SegmentationView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class SegmentationComponent extends LayoutComponent {
  static readonly typeName = "segmentation";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(SegmentationView), store, parent);
  }
}
