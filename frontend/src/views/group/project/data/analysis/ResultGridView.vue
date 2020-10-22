<template>
  <grid-layout
    :layout.sync="layout"
    :col-num="12"
    :row-height="600"
    :is-draggable="true"
    :is-resizable="true"
    :is-mirrored="false"
    :vertical-compact="true"
    :margin="[10, 10]"
    :use-css-transforms="true"
  >
    <grid-item
      v-for="item in layout"
      :x="item.x"
      :y="item.y"
      :w="item.w"
      :h="item.h"
      :i="item.i"
      :key="item.i"
      drag-allow-from=".widget-title"
    >
      <ScatterWidget v-if="item.i === `scatter`" />
      <PcaWidget v-else-if="item.i === `pca`" />
      <TsneWidget v-else-if="item.i === `tsne`" />
      <UmapWidget v-else-if="item.i === `umap`" />
    </grid-item>
  </grid-layout>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import VueGridLayout from "vue-grid-layout";
import PcaWidget from "@/views/group/project/data/analysis/pca/PcaWidget.vue";
import ScatterWidget from "@/views/group/project/data/analysis/scatter/ScatterWidget.vue";
import TsneWidget from "@/views/group/project/data/analysis/tsne/TsneWidget.vue";
import UmapWidget from "@/views/group/project/data/analysis/umap/UmapWidget.vue";

@Component({
  components: {
    UmapWidget,
    TsneWidget,
    ScatterWidget,
    PcaWidget,
    GridLayout: VueGridLayout.GridLayout,
    GridItem: VueGridLayout.GridItem,
  },
})
export default class ResultGridView extends Vue {
  readonly layout = [
    { x: 0, y: 0, w: 6, h: 1, i: "scatter" },
    { x: 6, y: 0, w: 6, h: 1, i: "pca" },
    { x: 0, y: 1, w: 6, h: 1, i: "tsne" },
    { x: 6, y: 1, w: 6, h: 1, i: "umap" },
  ];
}
</script>

<style scoped>
.vue-grid-layout {
  background: #eee;
}
.vue-grid-item:not(.vue-grid-placeholder) {
  background: #ccc;
}
.vue-grid-item .resizing {
  opacity: 0.9;
}
.vue-grid-item .static {
  background: #cce;
}
.vue-grid-item .text {
  font-size: 24px;
  text-align: center;
  position: absolute;
  top: 0;
  bottom: 0;
  left: 0;
  right: 0;
  margin: auto;
  height: 100%;
  width: 100%;
}
.vue-grid-item .no-drag {
  height: 100%;
  width: 100%;
}
.vue-grid-item .minMax {
  font-size: 12px;
}
.vue-grid-item .add {
  cursor: pointer;
}
.vue-draggable-handle {
  position: absolute;
  width: 20px;
  height: 20px;
  top: 0;
  left: 0;
  background: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='10' height='10'><circle cx='5' cy='5' r='5' fill='#999999'/></svg>")
    no-repeat;
  background-position: bottom right;
  padding: 0 8px 8px 0;
  background-repeat: no-repeat;
  background-origin: content-box;
  box-sizing: border-box;
  cursor: pointer;
}
</style>

<style>
.widget-title {
  margin: 0;
  padding: 0 0 0 10px;
}
</style>
