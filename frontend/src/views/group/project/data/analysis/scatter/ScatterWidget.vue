<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">Scatter</span>
    </v-card-title>
    <v-divider />
    <v-toolbar dense flat>
      x:
      <v-select
        :items="markers"
        v-model="markerX"
        label="Marker X"
        :rules="[required]"
        hide-details
        dense
        solo
        flat
        class="select-input"
      />
      y:
      <v-select
        :items="markers"
        v-model="markerY"
        label="Marker Y"
        :rules="[required]"
        hide-details
        dense
        solo
        flat
        class="select-input"
      />
      <!--      <v-switch v-model="showRegression" label="Show regression" hide-details inset dense />-->
    </v-toolbar>
    <ScatterPlot2d
      v-if="plotData"
      :ignore-selection="true"
      plot-id="scatterPlot"
      :data="plotData"
      title="Scatter Plot"
      :x-axis-title="markerX"
      :y-axis-title="markerY"
      class="plot"
    />
  </v-card>
</template>

<script lang="ts">
import { datasetsModule } from "@/modules/datasets";
import { required } from "@/utils/validators";
import { Component, Vue, Watch } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { resultsModule } from "@/modules/results";

@Component({
  components: { ScatterPlot2d },
})
export default class ScatterWidget extends Vue {
  readonly datasetContext = datasetsModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);

  readonly required = required;

  showRegression = false;
  markerX: string | null = null;
  markerY: string | null = null;

  get plotData() {
    return this.resultsContext.getters.scatterPlotData;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get markers() {
    return this.activeDataset && this.activeDataset.meta["channel_map"]
      ? Object.keys(this.activeDataset.meta["channel_map"])
      : [];
  }

  @Watch("markerX")
  markerXChanged(value) {
    this.submit();
  }

  @Watch("markerY")
  markerYChanged(value) {
    this.submit();
  }

  @Watch("showRegression")
  showRegressionChanged(value) {
    this.submit();
  }

  async submit() {
    if (this.markerX && this.markerY) {
      await this.resultsContext.actions.getScatterPlotData({
        markerX: this.markerX,
        markerY: this.markerY,
      });
    }
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
  height: calc(100% - 86px);
}
.select-input {
  max-width: 200px;
  margin-right: 20px;
}
</style>
