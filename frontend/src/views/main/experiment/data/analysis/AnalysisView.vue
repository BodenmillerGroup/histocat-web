<template>
  <v-tabs v-model="tab">
    <v-tab>Scatter Plot</v-tab>
    <v-tab>Box Plot</v-tab>
    <v-tab>PCA</v-tab>
    <v-tab>t-SNE</v-tab>
    <v-tab>UMAP</v-tab>
    <v-tab>PhenoGraph</v-tab>
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
      <TSNETab />
    </v-tab-item>
    <v-tab-item>
      <UMAPTab />
    </v-tab-item>
    <v-tab-item>
      <PhenoGraphTab />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import BoxPlotTab from "@/views/main/experiment/data/analysis/box/BoxPlotTab.vue";
import PCATab from "@/views/main/experiment/data/analysis/pca/PCATab.vue";
import PhenoGraphTab from "@/views/main/experiment/data/analysis/phenograph/PhenoGraphTab.vue";
import TSNETab from "@/views/main/experiment/data/analysis/tsne/TSNETab.vue";
import ScatterPlotTab from "@/views/main/experiment/data/analysis/scatter/ScatterPlotTab.vue";
import UMAPTab from "@/views/main/experiment/data/analysis/umap/UMAPTab.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    PhenoGraphTab,
    UMAPTab,
    PCATab,
    BoxPlotTab,
    ScatterPlotTab,
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
