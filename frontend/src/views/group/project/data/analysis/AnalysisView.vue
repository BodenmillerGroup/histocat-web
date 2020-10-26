<template>
  <v-tabs v-model="tab">
    <v-tab>Pipeline</v-tab>
    <v-tab>Result</v-tab>
    <v-tab>Plots</v-tab>
    <v-tab>Box Plot</v-tab>
    <v-tab>PhenoGraph</v-tab>
    <v-tab-item>
      <PipelineView />
    </v-tab-item>
    <v-tab-item>
      <ResultGridView />
    </v-tab-item>
    <v-tab-item>
      <PlotsGridView />
    </v-tab-item>
    <v-tab-item>
      <BoxPlotTab />
    </v-tab-item>
    <v-tab-item>
      <PhenoGraphTab />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import BoxPlotTab from "@/views/group/project/data/analysis/box/BoxPlotTab.vue";
import PhenoGraphTab from "@/views/group/project/data/analysis/phenograph/PhenoGraphTab.vue";
import { Component, Vue } from "vue-property-decorator";
import PipelineView from "@/views/group/project/data/analysis/pipeline/PipelineView.vue";
import ResultGridView from "@/views/group/project/data/analysis/ResultGridView.vue";
import PlotsGridView from "@/views/group/project/data/analysis/PlotsGridView.vue";

@Component({
  components: {
    PlotsGridView,
    ResultGridView,
    PipelineView,
    PhenoGraphTab,
    BoxPlotTab,
  },
})
export default class AnalysisView extends Vue {
  tab = 0;
}
</script>
