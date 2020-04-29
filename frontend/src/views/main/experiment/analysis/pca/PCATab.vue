<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <Scatter2D v-if="nComponents === '2'" :data="pcaData" title="Principal Component Analysis" :width="responsive.width - 500" :height="responsive.height - 150" />
      <Scatter3D v-else :data="pcaData" title="Principal Component Analysis" />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>PCA Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
            <v-chip v-for="item in channels" :key="item" :value="item" small>
              {{ item }}
            </v-chip>
          </v-chip-group>
          <v-card-actions>
            <v-btn @click="selectAll" small :disabled="selectedChannels.length === channels.length">
              Select all
            </v-btn>
            <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0">
              Clear all
            </v-btn>
          </v-card-actions>
          <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
            <v-radio label="2D" value="2" />
            <v-radio label="3D" value="3" />
          </v-radio-group>
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
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import Scatter3D from "@/components/charts/Scatter3D.vue";
import Scatter2D from "@/components/charts/Scatter2D.vue";
import {responsiveModule} from "@/modules/responsive";

@Component({
  components: { Scatter2D, Scatter3D },
})
export default class PCATab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly responsiveContext = responsiveModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  selectedChannels: any[] = [];
  nComponents = "2";
  heatmap: { type: string; label: string } | null = null;

  points: any[] = [];

  get responsive() {
    return this.responsiveContext.getters.responsive;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    return this.showOptions ? 9 : 12;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get channels() {
    return this.datasetContext.getters.channels;
  }

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  selectAll() {
    this.selectedChannels = this.channels;
  }

  clearAll() {
    this.selectedChannels = [];
  }

  async submit() {
    let heatmap = "";
    if (this.heatmap) {
      heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
    }

    const acquisitionIds =
      this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

    await this.analysisContext.actions.getPCAData({
      dataset_id: this.activeDataset!.id,
      acquisition_ids: acquisitionIds,
      n_components: parseInt(this.nComponents, 10),
      heatmapType: this.heatmap ? this.heatmap.type : "",
      heatmap: heatmap,
      markers: this.selectedChannels,
    });
  }

  get pcaData() {
    return this.analysisContext.getters.pcaData;
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
