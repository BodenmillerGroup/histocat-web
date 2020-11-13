<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">PCA</span>
    </v-card-title>
    <v-divider />
    <ScatterPlot2d v-if="plotData" plot-id="pcaPlot" :ignore-selection="false" :data="plotData" title="PCA" class="plot" />
  </v-card>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { resultsModule } from "@/modules/results";

@Component({
  components: { ScatterPlot2d },
})
export default class PcaWidget extends Vue {
  readonly resultsContext = resultsModule.context(this.$store);

  get plotData() {
    return this.resultsContext.getters.pcaData;
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
