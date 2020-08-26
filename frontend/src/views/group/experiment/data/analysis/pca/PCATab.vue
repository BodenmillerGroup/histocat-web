<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">Please select dataset</v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <div v-else :class="layoutClass">
    <Scatter2D v-if="nComponents === '2'" plot-id="pca2D" :data="pcaData" title="2D PCA" />
    <Scatter3D v-else plot-id="pca3D" :data="pcaData" :show-regression="false" title="3D PCA" />
    <div v-if="showOptions">
      <v-card tile>
        <v-card-title>PCA Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
            <v-chip v-for="channel in channels" :key="channel" :value="channel" small>
              {{ channel }}
            </v-chip>
          </v-chip-group>
          <v-card-actions>
            <v-btn @click="selectAll" small :disabled="selectedChannels.length === channels.length"> Select all </v-btn>
            <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0"> Clear all </v-btn>
          </v-card-actions>
          <v-row>
            <v-col>
              <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
                <v-radio label="2D" value="2" />
                <v-radio label="3D" value="3" />
              </v-radio-group>
            </v-col>
          </v-row>
          <v-row>
            <v-col>
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
            </v-col>
          </v-row>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0"> Analyze </v-btn>
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
import { Component, Vue } from "vue-property-decorator";
import Scatter3D from "@/components/charts/Scatter3D.vue";
import Scatter2D from "@/components/charts/Scatter2D.vue";

@Component({
  components: { Scatter2D, Scatter3D },
})
export default class PCATab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  selectedChannels: any[] = [];
  nComponents = "2";
  heatmap: { type: string; label: string } | null = null;

  get pcaData() {
    return this.analysisContext.getters.pcaData;
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
    const heatmap = this.heatmap ? this.heatmap.label : "";
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
