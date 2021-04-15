<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">PCA</span>
    </v-card-title>
    <v-divider />
    <ScatterPlot2d
      v-if="hasData"
      plot-id="pcaPlot"
      :ignore-selection="false"
      :data="plotData"
      mapping="pca"
      title="PCA"
      x-axis-title="PC1"
      y-axis-title="PC2"
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
export default class PcaWidget extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get hasData() {
    return this.cellsContext.getters.activeResult && this.cellsContext.getters.activeResult.output.pca;
  }

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
