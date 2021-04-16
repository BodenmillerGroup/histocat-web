<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">UMAP</span>
    </v-card-title>
    <v-divider />
    <ScatterPlot2d
      plot-id="umapPlot"
      :ignore-selection="false"
      :data="plotData"
      mapping="umap"
      title="UMAP"
      x-axis-title="UMAP1"
      y-axis-title="UMAP2"
      class="plot"
    />
  </v-card>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { cellsModule } from "@/modules/cells";

@Component({
  components: { ScatterPlot2d },
})
export default class UmapWidget extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get plotData() {
    return this.cellsContext.getters.cellsByAcquisition;
  }
}
</script>

<style scoped>
.widget-container {
  width: 100%;
  height: 100%;
}
.plot {
  width: 100%;
  height: calc(100% - 36px);
}
</style>
