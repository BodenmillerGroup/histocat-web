import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import LouvainView from "./LouvainView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class LouvainComponent extends LayoutComponent {
  static readonly typeName = "louvain";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(LouvainView), store, parent);
  }
}
