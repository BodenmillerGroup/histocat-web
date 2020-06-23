<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisitionId" icon="mdi-alert-circle-outline">
    Please select acquisition
  </v-banner>
  <div v-else :class="layoutClass">
    <BoxPlotView plot-id="boxPlot" :data="boxPlotData" title="BoxPlot" />
    <div v-if="showOptions">
      <v-card tile>
        <v-card-title>Box Plot Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedItems" multiple column active-class="primary--text">
            <v-chip v-for="item in items" :key="item" :value="item" small>
              {{ item }}
            </v-chip>
          </v-chip-group>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="selectAll" small :disabled="selectedItems.length === items.length">
            Select all
          </v-btn>
          <v-btn @click="clearAll" small :disabled="selectedItems.length === 0">
            Clear all
          </v-btn>
        </v-card-actions>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedItems.length === 0">
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
import { Component, Vue } from "vue-property-decorator";
import BoxPlotView from "@/views/main/experiment/data/analysis/box/BoxPlotView.vue";

@Component({
  components: { BoxPlotView },
})
export default class BoxPlotTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  selectedItems: any[] = [];

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

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get items() {
    return this.activeDataset && this.activeDataset.input["channel_map"]
      ? Object.keys(this.activeDataset.input["channel_map"])
      : [];
  }

  selectAll() {
    this.selectedItems = this.items;
  }

  clearAll() {
    this.selectedItems = [];
  }

  async submit() {
    if (!this.activeDataset) {
      self.alert("Please select a dataset");
      return;
    }

    if (!this.activeAcquisitionId) {
      self.alert("Please select an acquisition");
      return;
    }

    await this.analysisContext.actions.getBoxPlotData({
      datasetId: this.activeDataset.id,
      acquisitionId: this.activeAcquisitionId,
      markers: this.selectedItems,
    });
  }

  get boxPlotData() {
    return this.analysisContext.getters.boxPlotData;
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
