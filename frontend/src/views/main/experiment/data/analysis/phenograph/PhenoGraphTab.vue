<template>
  <v-banner v-if="!activeDatasetId" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisitionId && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <div v-else :class="layoutClass">
    <PhenoGraphView plot-id="phenographHeatmap" :data="phenographData" title="PhenoGraph Clustering" />
    <div v-if="showOptions">
      <v-card tile>
        <v-card-title>PhenoGraph Settings</v-card-title>
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
              <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0">
                Clear all
              </v-btn>
            </v-card-actions>
            <v-radio-group v-model="jaccard" mandatory hide-details label="Mode">
              <v-radio label="Jaccard metric" value="jaccard" />
              <v-radio label="Gaussian kernel" value="gaussian" />
            </v-radio-group>
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  max="200"
                  step="1"
                  label="Neighbors"
                  v-model.number="nearestNeighbors"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  step="1"
                  label="Minimum cluster size"
                  v-model.number="minClusterSize"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
            <v-select :items="metrics" v-model="primaryMetric" label="Metric" hide-details dense class="mt-5" />
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
        <v-card-text>
          <v-select
            :items="results"
            v-model="result"
            item-text="name"
            return-object
            label="Results"
            hint="PhenoGraph processed data"
            persistent-hint
            clearable
            dense
            @change="resultChanged"
          />
        </v-card-text>
        <v-card-actions>
          <v-btn @click="display" color="primary" block :disabled="!result">
            Display
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
import { required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SELECTED_ACQUISITION_IDS } from "@/modules/experiment/events";
import PhenoGraphView from "@/views/main/experiment/data/analysis/phenograph/PhenoGraphView.vue";

@Component({
  components: { PhenoGraphView },
})
export default class PhenoGraphTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  readonly required = required;
  readonly metrics = ["euclidean", "manhattan", "correlation", "cosine"];

  valid = false;

  selectedChannels: any[] = [];
  nearestNeighbors = 30;
  jaccard = "jaccard";
  minClusterSize = 10;
  primaryMetric = "euclidean";

  result: any = null;

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get layoutClass() {
    if (!this.showOptions) {
      return "layout-without-options";
    }
    return "layout-full";
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  get activeDatasetId() {
    return this.datasetContext.getters.activeDatasetId;
  }

  get channels() {
    return this.datasetContext.getters.channels;
  }

  get results() {
    return this.datasetContext.getters.phenographResults;
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
        this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisitionId!];

      await this.analysisContext.actions.submitPhenoGraph({
        dataset_id: this.activeDatasetId!,
        acquisition_ids: acquisitionIds,
        markers: this.selectedChannels,
        jaccard: this.jaccard === "jaccard",
        min_cluster_size: this.minClusterSize,
        nearest_neighbors: this.nearestNeighbors,
        primary_metric: this.primaryMetric,
      });
    }
  }

  resultChanged(result) {
    if (result) {
      this.selectedChannels = result.params.markers;
      this.nearestNeighbors = result.params.nearest_neighbors;
      this.primaryMetric = result.params.primary_metric;
      this.minClusterSize = result.params.min_cluster_size;
      this.jaccard = result.params.jaccard ? "jaccard" : "gaussian";

      BroadcastManager.publish(SET_SELECTED_ACQUISITION_IDS, result.params.acquisition_ids);
    }
  }

  async display() {
    await this.analysisContext.actions.getPhenoGraphResult({
      datasetId: this.activeDatasetId!,
      name: this.result ? this.result.name : "",
    });
  }

  get phenographData() {
    return this.analysisContext.getters.phenographData;
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
