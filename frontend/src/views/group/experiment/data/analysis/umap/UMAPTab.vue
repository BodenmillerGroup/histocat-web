<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">Please select dataset</v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <div v-else :class="layoutClass">
    <Scatter2D v-if="nComponents === '2'" plot-id="umap2D" :data="umapData" title="2D UMAP" />
    <Scatter3D v-else plot-id="umap3D" :data="umapData" title="3D UMAP" />
    <div v-if="showOptions">
      <v-card tile>
        <v-card-title>UMAP Settings</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
            <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
              <v-chip v-for="item in channels" :key="item" :value="item" small>
                {{ item }}
              </v-chip>
            </v-chip-group>
            <v-card-actions>
              <v-btn @click="selectAll" small :disabled="selectedChannels.length === channels.length">
                Select all
              </v-btn>
              <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0"> Clear all </v-btn>
            </v-card-actions>
            <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
              <v-radio label="2D" value="2" />
              <v-radio label="3D" value="3" />
            </v-radio-group>
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  max="200"
                  step="1"
                  label="Neighbors"
                  v-model.number="nNeighbors"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="0.0"
                  max="0.99"
                  step="0.01"
                  label="Minimum distance"
                  v-model.number="minDist"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
            <v-select :items="metrics" v-model="metric" label="Metric" dense hide-details class="mt-5" />
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0"> Analyze </v-btn>
        </v-card-actions>
        <v-card-text>
          <v-select
            :items="results"
            v-model="result"
            item-text="name"
            return-object
            label="Results"
            hint="UMAP processed data"
            persistent-hint
            clearable
            dense
            @change="resultChanged"
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
        </v-card-text>
        <v-card-actions>
          <v-btn @click="display" color="primary" block :disabled="!result"> Display </v-btn>
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
import Scatter2D from "@/components/charts/Scatter2D.vue";
import Scatter3D from "@/components/charts/Scatter3D.vue";
import { resultModule } from "@/modules/results";

@Component({
  components: { Scatter3D, Scatter2D },
})
export default class UMAPTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly resultContext = resultModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;
  readonly metrics = [
    "euclidean",
    "manhattan",
    "chebyshev",
    "minkowski",
    "canberra",
    "braycurtis",
    "haversine",
    "mahalanobis",
    "wminkowski",
    "seuclidean",
    "cosine",
    "correlation",
  ];

  valid = false;

  selectedChannels: any[] = [];
  nComponents = "2";
  nNeighbors = 15;
  minDist = 0.1;
  metric = "euclidean";

  heatmap: { type: string; label: string } | null = null;

  result: any = null;

  get umapData() {
    return this.analysisContext.getters.umapData;
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

  get results() {
    return this.resultContext.getters.umapResults;
  }

  selectAll() {
    this.selectedChannels = this.channels;
  }

  clearAll() {
    this.selectedChannels = [];
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const acquisitionIds =
        this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

      await this.analysisContext.actions.submitUMAP({
        dataset_id: this.activeDataset!.id,
        acquisition_ids: acquisitionIds,
        n_components: parseInt(this.nComponents, 10),
        markers: this.selectedChannels,
        n_neighbors: this.nNeighbors,
        min_dist: this.minDist,
        metric: this.metric,
      });
    }
  }

  resultChanged(result) {
    if (result) {
      this.nComponents = result.params.n_components.toString();
      this.nNeighbors = result.params.n_neighbors;
      this.minDist = result.params.min_dist;
      this.metric = result.params.metric;
      this.selectedChannels = result.params.markers;

      this.experimentContext.actions.setSelectedAcquisitionIds(result.params.acquisition_ids);
    }
  }

  async display() {
    if (!this.activeDataset) {
      self.alert("Please select a dataset");
      return;
    }

    if (!this.result) {
      self.alert("Please select result data");
      return;
    }

    const heatmap = this.heatmap ? this.heatmap.label : "";

    if (this.result) {
      await this.analysisContext.actions.getUMAPResult({
        resultId: this.result.id,
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
