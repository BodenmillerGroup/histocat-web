<template>
  <v-tabs v-model="tab">
    <v-tab>Pipeline</v-tab>
    <v-tab>Result</v-tab>
    <v-tab v-if="activeResult && (activeResult.output.leiden || activeResult.output.louvain)">Clustering</v-tab>
    <v-tab-item>
      <PipelineView />
    </v-tab-item>
    <v-tab-item>
      <ResultGridView />
    </v-tab-item>
    <v-tab-item v-if="activeResult && (activeResult.output.leiden || activeResult.output.louvain)">
      <ClusteringGridView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import PipelineView from "@/views/group/project/data/analysis/pipeline/PipelineView.vue";
import ResultGridView from "@/views/group/project/data/analysis/ResultGridView.vue";
import { resultsModule } from "@/modules/results";
import ClusteringGridView from "@/views/group/project/data/analysis/clustering/ClusteringGridView.vue";

@Component({
  components: {
    ClusteringGridView,
    ResultGridView,
    PipelineView,
  },
})
export default class AnalysisView extends Vue {
  readonly resultsContext = resultsModule.context(this.$store);

  tab = 0;

  get activeResult() {
    return this.resultsContext.getters.activeResult;
  }
}
</script>
