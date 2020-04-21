<template>
  <v-row no-gutters class="chart-container">
    <v-col :cols="columns">
      <Graph :key="graphRenderCounter" />
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Component, Vue, Watch } from "vue-property-decorator";
import Graph from "@/cellxgene/components/graph/graph.vue";
import {controlsModule} from "@/modules/controls";

@Component({
  components: { Graph },
})
export default class TestTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly controlsContext = controlsModule.context(this.$store);

  get graphRenderCounter() {
    return this.controlsContext.getters.graphRenderCounter;
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

  async submit() {
    const acquisitionIds =
      this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

    // await this.analysisContext.actions.getPCAData({
    //   dataset_id: this.activeDataset!.id,
    //   acquisition_ids: acquisitionIds,
    //   n_components: parseInt(this.nComponents, 10),
    //   heatmapType: this.heatmap ? this.heatmap.type : "",
    //   heatmap: heatmap,
    //   markers: this.selectedChannels,
    // });
  }

  @Watch("data")
  pcaDataChanged(data) {
    if (data) {
    }
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
