<template>
  <v-tabs v-model="tab">
    <v-tab>Pipeline</v-tab>
    <v-tab>Scatter Plot</v-tab>
    <v-tab>Box Plot</v-tab>
    <v-tab>PCA</v-tab>
    <v-tab>t-SNE</v-tab>
    <v-tab>UMAP</v-tab>
    <v-tab>PhenoGraph</v-tab>
    <v-tab-item>
      <PipelineView />
    </v-tab-item>
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
import BoxPlotTab from "@/views/group/project/data/analysis/box/BoxPlotTab.vue";
import PCATab from "@/views/group/project/data/analysis/pca/PCATab.vue";
import PhenoGraphTab from "@/views/group/project/data/analysis/phenograph/PhenoGraphTab.vue";
import TSNETab from "@/views/group/project/data/analysis/tsne/TSNETab.vue";
import ScatterPlotTab from "@/views/group/project/data/analysis/scatter/ScatterPlotTab.vue";
import UMAPTab from "@/views/group/project/data/analysis/umap/UMAPTab.vue";
import { Component, Vue } from "vue-property-decorator";
import PipelineView from "@/views/group/project/data/analysis/pipeline/PipelineView.vue";

@Component({
  components: {
    PipelineView,
    ScatterPlotTab,
    PhenoGraphTab,
    UMAPTab,
    PCATab,
    BoxPlotTab,
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
