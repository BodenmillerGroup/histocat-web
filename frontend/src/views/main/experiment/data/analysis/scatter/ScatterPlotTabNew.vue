<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <div v-else :class="layoutClass">
    <ScatterPlotView
      v-if="!markerZ"
      :data="scatterPlotData"
      :show-regression="showRegression"
      :regression-type="regressionType"
      :polynomial-order="polynomialOrder"
      title="2D Scatter Plot"
    />
    <Scatter3D v-else :data="scatterPlotData" title="3D Scatter Plot" />
    <div v-if="showOptions">
      <v-card tile>
        <v-card-title>Scatter Plot Settings</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
            <v-select
              :items="items"
              v-model="markerX"
              label="X"
              hint="X axis marker"
              persistent-hint
              :rules="[required]"
              dense
            />
            <v-select
              :items="items"
              v-model="markerY"
              label="Y"
              hint="Y axis marker"
              persistent-hint
              :rules="[required]"
              dense
              class="mt-5"
            />
            <v-select
              :items="items"
              v-model="markerZ"
              label="Z"
              hint="Z axis marker"
              persistent-hint
              clearable
              dense
              class="mt-5"
            />
            <v-select
              :items="heatmaps"
              v-model="heatmap"
              label="Heatmap"
              hint="Heatmap marker"
              item-text="label"
              return-object
              persistent-hint
              clearable
              dense
              class="mt-5"
            />
            <v-switch v-if="!markerZ" v-model="showRegression" label="Show regression" />
            <v-select
              v-if="!markerZ"
              :items="regressionTypes"
              v-model="regressionType"
              label="Regression type"
              hide-details
              dense
              class="mt-5"
            />
            <v-text-field
              v-if="!markerZ && regressionType === 'polynomial'"
              type="number"
              min="2"
              step="1"
              label="Polynomial order"
              v-model.number="polynomialOrder"
              :rules="[required]"
              hide-details
              class="mt-5"
            />
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="!valid">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </div>
  </div>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";
import ScatterPlotView from "@/views/main/experiment/data/analysis/scatter/ScatterPlotView.vue";
import Scatter3D from "@/components/charts/Scatter3D.vue";

type RegressionType = "linear" | "polynomial";

@Component({
  components: { Scatter3D, ScatterPlotView },
})
export default class ScatterPlotTabNew extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;
  readonly regressionTypes: RegressionType[] = ["linear", "polynomial"];

  valid = false;

  showRegression = false;
  regressionType: RegressionType = "linear";
  polynomialOrder = 2;

  markerX: string | null = null;
  markerY: string | null = null;
  markerZ: string | null = null;
  heatmap: { type: string; label: string } | null = null;

  get scatterPlotData() {
    return this.analysisContext.getters.scatterPlotData;
  }

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get layoutClass() {
    if (!this.showOptions) {
      return "layout-without-options";
    }
    return "layout-full";
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  get items() {
    return this.activeDataset && this.activeDataset.input["channel_map"]
      ? Object.keys(this.activeDataset.input["channel_map"])
      : [];
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      let heatmap = "";
      if (this.heatmap) {
        heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
      }

      const acquisitionIds =
        this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

      await this.analysisContext.actions.getScatterPlotData({
        datasetId: this.activeDataset!.id,
        acquisitionIds: acquisitionIds,
        markerX: this.markerX!,
        markerY: this.markerY!,
        markerZ: this.markerZ ? this.markerZ : "",
        heatmapType: this.heatmap ? this.heatmap.type : "",
        heatmap: heatmap,
      });
    }
  }
}
</script>

<style scoped>
.layout-full {
  display: grid;
  grid-template-columns: 1fr 380px;
  grid-template-rows: auto;
}
.layout-without-options {
  display: grid;
  grid-template-columns: 1fr;
  grid-template-rows: auto;
}
</style>
