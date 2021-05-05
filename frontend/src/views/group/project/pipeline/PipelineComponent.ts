import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import PipelineView from "./PipelineView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class PipelineComponent extends LayoutComponent {
  static readonly typeName = PipelineComponent.name;

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(PipelineView), store, parent);
  }
}
