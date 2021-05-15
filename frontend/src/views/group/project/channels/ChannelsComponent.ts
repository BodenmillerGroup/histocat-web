import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import ChannelsView from "./ChannelsView.vue";
import { LayoutComponent } from "@/components/LayoutComponent";

export class ChannelsComponent extends LayoutComponent {
  static readonly typeName = "channels";

  constructor(_container: ComponentContainer, store: any, parent) {
    super(_container, Vue.extend(ChannelsView), store, parent);
  }
}
