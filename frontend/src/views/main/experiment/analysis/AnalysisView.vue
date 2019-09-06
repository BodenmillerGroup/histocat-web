<template>
  <v-row no-gutters>
    <v-col>
      <v-tabs v-model="tab">
        <v-tab>Scatter Plot</v-tab>
        <v-tab>Box Plot</v-tab>
        <v-tab>PCA</v-tab>
        <v-tab>t-SNE</v-tab>
        <v-tab>Segmentation</v-tab>
        <v-tab-item>
          <ScatterPlotTab/>
        </v-tab-item>
        <v-tab-item>
          <BoxPlotTab/>
        </v-tab-item>
        <v-tab-item>
          <PCATab/>
        </v-tab-item>
        <v-tab-item>
          <TSNETab/>
        </v-tab-item>
        <v-tab-item>
          <SegmentationTab/>
        </v-tab-item>
      </v-tabs>
    </v-col>
  </v-row>
</template>

<script lang="ts">
  import BoxPlotTab from '@/views/main/experiment/analysis/box/BoxPlotTab.vue';
  import PCATab from '@/views/main/experiment/analysis/pca/PCATab.vue';
  import TSNETab from '@/views/main/experiment/analysis/tsne/TSNETab.vue';
  import ScatterPlotTab from '@/views/main/experiment/analysis/scatter/ScatterPlotTab.vue';
  import SegmentationTab from '@/views/main/experiment/analysis/segmentation/SegmentationTab.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
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
