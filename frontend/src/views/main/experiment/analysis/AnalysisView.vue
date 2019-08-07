<template>
  <v-layout row>
    <v-flex>
      <v-tabs v-model="tab">
        <v-tab>Segmentation</v-tab>
        <v-tab-item>
          <SegmentationTab/>
        </v-tab-item>
      </v-tabs>
    </v-flex>
  </v-layout>
</template>

<script lang="ts">
  import SegmentationTab from '@/views/main/experiment/analysis/segmentation/SegmentationTab.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      SegmentationTab,
    },
  })
  export default class AnalysisView extends Vue {
    tab = 0;
  }
</script>
