<template>
  <div class="widget-container">
    <v-toolbar v-if="activeResult" dense flat>
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
        item-text="label"
        item-value="tag"
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
        item-text="label"
        item-value="tag"
      />
    </v-toolbar>
    <ScatterPlot2d
      v-if="activeResult"
      :ignore-selection="true"
      plot-id="scatterPlot"
      :data="plotData"
      mapping="scatterplot"
      title="Scatter Plot"
      :x-axis-title="markerX"
      :y-axis-title="markerY"
      class="plot"
    />
  </div>
</template>

<script lang="ts">
import { required } from "@/utils/validators";
import { Component, Vue, Watch } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { projectsModule } from "@/modules/projects";
import { cellsModule } from "@/modules/cells";

@Component({
  components: { ScatterPlot2d },
})
export default class ScatterView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  readonly required = required;

  markerX: string | null = null;
  markerY: string | null = null;

  get activeResult() {
    return this.cellsContext.getters.activeResult;
  }

  get plotData() {
    return this.cellsContext.getters.cellsByAcquisition;
  }

  get channels() {
    const acquisition = this.projectsContext.getters.activeAcquisition;
    return acquisition ? acquisition.channels : null;
  }

  get markers() {
    return this.cellsContext.getters.markers.map((tag) => {
      return {
        tag: tag,
        label: this.channels && this.channels[tag] ? this.channels[tag].customLabel : tag,
      };
    });
  }

  @Watch("markers")
  markersChanged(value) {
    this.markerX = null;
    this.markerY = null;
  }

  @Watch("markerX")
  markerXChanged(value) {
    this.submit();
  }

  @Watch("markerY")
  markerYChanged(value) {
    this.submit();
  }

  async submit() {
    if (this.markerX && this.markerY) {
      await this.cellsContext.actions.getScatterPlotData({
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
