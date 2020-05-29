<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <Scatter2D
        v-if="nComponents === '2'"
        :data="tsneData"
        title="t-Distributed Stochastic Neighbor Embedding"
        :width="responsive.width - 500"
        :height="responsive.height - 150"
      />
      <Scatter3D v-else :data="tsneData" title="t-Distributed Stochastic Neighbor Embedding" />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>t-SNE Settings</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
            <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
              <v-chip v-for="channel in channels" :key="channel" :value="channel" small>
                {{ channel }}
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
            <v-row>
              <v-col>
                <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
                  <v-radio label="2D" value="2" />
                  <v-radio label="3D" value="3" />
                </v-radio-group>
              </v-col>
              <v-col>
                <v-radio-group v-model="init" mandatory hide-details label="Initialization">
                  <v-radio label="PCA" value="pca" />
                  <v-radio label="Random" value="random" />
                </v-radio-group>
              </v-col>
            </v-row>
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="5"
                  max="50"
                  step="1"
                  label="Perplexity"
                  v-model.number="perplexity"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="10"
                  max="1000"
                  step="1"
                  label="Learning rate"
                  v-model.number="learningRate"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="250"
                  step="1"
                  label="Iterations"
                  v-model.number="iterations"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="0"
                  max="1"
                  step="0.01"
                  label="Angle (theta)"
                  v-model.number="theta"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
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
            hint="t-SNE processed data"
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
          <v-btn @click="display" color="primary" block :disabled="!result">
            Display
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
import { required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";
import Scatter2D from "@/components/charts/Scatter2D.vue";
import Scatter3D from "@/components/charts/Scatter3D.vue";
import { responsiveModule } from "@/modules/responsive";
import {BroadcastManager} from "@/utils/BroadcastManager";
import {SET_SELECTED_ACQUISITION_IDS} from "@/modules/experiment/events";

@Component({
  components: { Scatter3D, Scatter2D },
})
export default class TSNETab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly responsiveContext = responsiveModule.context(this.$store);

  readonly required = required;

  valid = false;

  selectedChannels: any[] = [];
  nComponents = "2";
  perplexity = 30;
  learningRate = 200;
  iterations = 1000;
  theta = 0.5;
  init = "pca";

  heatmap: { type: string; label: string } | null = null;

  result: any = null;

  get tsneData() {
    return this.analysisContext.getters.tsneData;
  }

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

  get results() {
    return this.datasetContext.getters.tsneResults;
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

      await this.analysisContext.actions.submitTSNE({
        dataset_id: this.activeDataset!.id,
        acquisition_ids: acquisitionIds,
        n_components: parseInt(this.nComponents, 10),
        markers: this.selectedChannels,
        perplexity: this.perplexity,
        learning_rate: this.learningRate,
        iterations: this.iterations,
        theta: this.theta,
        init: this.init,
      });
    }
  }

  resultChanged(result) {
    if (result) {
      this.nComponents = result.params.n_components.toString();
      this.perplexity = result.params.perplexity;
      this.iterations = result.params.iterations;
      this.learningRate = result.params.learning_rate;
      this.selectedChannels = result.params.markers;
      this.theta = result.params.theta;
      this.init = result.params.init;

      BroadcastManager.publish(SET_SELECTED_ACQUISITION_IDS, result.params.acquisition_ids);
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

    let heatmap = "";
    if (this.heatmap) {
      heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
    }

    await this.analysisContext.actions.getTSNEResult({
      datasetId: this.activeDataset.id,
      name: this.result ? this.result.name : "",
      heatmapType: this.heatmap ? this.heatmap.type : "",
      heatmap: heatmap,
    });
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
