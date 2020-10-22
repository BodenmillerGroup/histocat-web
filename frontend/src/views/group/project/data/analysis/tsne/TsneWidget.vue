<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">tSNE</span>
    </v-card-title>
    <v-divider />
    <ScatterPlot2d plot-id="tsnePlot" :data="plotData" title="tSNE" class="plot" />
  </v-card>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { analysisModule } from "@/modules/analysis";
import { resultsModule } from "@/modules/results";

@Component({
  components: { ScatterPlot2d },
})
export default class TsneWidget extends Vue {
  readonly analysisContext = analysisModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);

  get heatmap() {
    return this.resultsContext.getters.heatmap;
  }

  get plotData() {
    return this.analysisContext.getters.tsneData;
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
