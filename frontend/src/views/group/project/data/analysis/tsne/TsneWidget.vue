<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">tSNE</span>
    </v-card-title>
    <v-divider />
    <ScatterPlot2d
      v-if="hasData"
      plot-id="tsnePlot"
      :ignore-selection="false"
      :data="plotData"
      mapping="tsne"
      title="tSNE"
      x-axis-title="tSNE1"
      y-axis-title="tSNE2"
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
export default class TsneWidget extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get hasData() {
    return this.cellsContext.getters.activeResult && this.cellsContext.getters.activeResult.output.tsne;
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
