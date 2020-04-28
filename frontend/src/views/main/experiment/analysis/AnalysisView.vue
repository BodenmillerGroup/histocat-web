<template>
  <v-row no-gutters>
    <v-col>
      <v-tabs v-model="tab">
        <v-tab>Scatter Plot</v-tab>
        <v-tab>Box Plot</v-tab>
        <v-tab>PCA</v-tab>
        <v-tab>PCA2</v-tab>
        <v-tab>t-SNE</v-tab>
        <v-tab>UMAP</v-tab>
        <v-tab>PhenoGraph</v-tab>
        <v-tab>Segmentation</v-tab>
        <v-tab>Test</v-tab>
        <v-tab-item>
          <ScatterPlotTab />
        </v-tab-item>
        <v-tab-item>
          <BoxPlotTab />
        </v-tab-item>
        <v-tab-item>
          <PCATab />
        </v-tab-item>
        <v-tab-item>
          <PCATab2 />
        </v-tab-item>
        <v-tab-item>
          <TSNETab />
        </v-tab-item>
        <v-tab-item>
          <UMAPTab />
        </v-tab-item>
        <v-tab-item>
          <PhenoGraphTab />
        </v-tab-item>
        <v-tab-item>
          <SegmentationTab />
        </v-tab-item>
        <v-tab-item>
          <TestTab />
        </v-tab-item>
      </v-tabs>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import BoxPlotTab from "@/views/main/experiment/analysis/box/BoxPlotTab.vue";
import PCATab from "@/views/main/experiment/analysis/pca/PCATab.vue";
import PhenoGraphTab from "@/views/main/experiment/analysis/phenograph/PhenoGraphTab.vue";
import TSNETab from "@/views/main/experiment/analysis/tsne/TSNETab.vue";
import ScatterPlotTab from "@/views/main/experiment/analysis/scatter/ScatterPlotTab.vue";
import SegmentationTab from "@/views/main/experiment/analysis/segmentation/SegmentationTab.vue";
import UMAPTab from "@/views/main/experiment/analysis/umap/UMAPTab.vue";
import { Component, Vue } from "vue-property-decorator";
import TestTab from "@/views/main/experiment/analysis/test/TestTab.vue";
import PCATab2 from "@/views/main/experiment/analysis/pca2/PCATab2.vue";

@Component({
  components: {
    PCATab2,
    TestTab,
    PhenoGraphTab,
    UMAPTab,
    PCATab,
    BoxPlotTab,
    ScatterPlotTab,
    SegmentationTab,
    TSNETab,
  },
})
export default class AnalysisView extends Vue {
  tab = 0;
}
</script>

<style>
/**
   * The default size is 600px×400px, for responsive charts
   * you may need to set percentage values as follows (also
   * don't forget to provide a size for the container).
   */
.echarts {
  width: 100%;
  height: 100%;
}

.v-slide-group__prev {
  min-width: 0;
}
</style>
