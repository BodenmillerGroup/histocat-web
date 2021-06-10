<template>
  <div class="widget-container">
    <ScatterPlot2d
      plot-id="umapPlot"
      ref="umapPlot"
      :ignore-selection="false"
      :data="cellsList"
      mapping="umap"
      title="UMAP"
      x-axis-title="UMAP1"
      y-axis-title="UMAP2"
    />
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { cellsModule } from "@/modules/cells";

@Component({
  components: { ScatterPlot2d },
})
export default class UmapView extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get cellsList() {
    return this.cellsContext.getters.cellsList;
  }

  refresh() {
    if (this.$refs.umapPlot) {
      (this.$refs.umapPlot as any).refresh();
    }
  }
}
</script>

<style scoped>
.widget-container {
  width: 100%;
  height: 100%;
}
</style>
