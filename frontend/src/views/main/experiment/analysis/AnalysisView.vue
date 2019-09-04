<template>
  <v-layout row>
    <v-flex>
      <v-tabs v-model="tab">
        <v-tab>Scatter Plot</v-tab>
        <v-tab>Box Plot</v-tab>
        <v-tab>Segmentation</v-tab>
        <v-tab-item>
          <ScatterPlotTab/>
        </v-tab-item>
        <v-tab-item>
          <BoxPlotTab/>
        </v-tab-item>
        <v-tab-item>
          <SegmentationTab/>
        </v-tab-item>
      </v-tabs>
    </v-flex>
  </v-layout>
</template>

<script lang="ts">
  import BoxPlotTab from '@/views/main/experiment/analysis/box/BoxPlotTab.vue';
  import ScatterPlotTab from '@/views/main/experiment/analysis/scatter/ScatterPlotTab.vue';
  import SegmentationTab from '@/views/main/experiment/analysis/segmentation/SegmentationTab.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      BoxPlotTab,
      ScatterPlotTab,
      SegmentationTab,
    },
  })
  export default class AnalysisView extends Vue {
    tab = 0;
  }
</script>
