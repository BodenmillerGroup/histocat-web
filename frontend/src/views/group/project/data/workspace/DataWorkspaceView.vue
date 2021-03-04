<template>
  <v-tabs v-model="tab">
    <v-tab>Datasets</v-tab>
    <v-tab>Results</v-tab>
    <v-tab>Pipelines</v-tab>
    <v-tab>Gates</v-tab>
    <v-tab-item>
      <DatasetsView />
    </v-tab-item>
    <v-tab-item>
      <ResultsView />
    </v-tab-item>
    <v-tab-item>
      <PipelinesView />
    </v-tab-item>
    <v-tab-item>
      <GatesView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import DatasetsView from "@/views/group/project/data/workspace/DatasetsView.vue";
import GatesView from "@/views/group/project/data/workspace/GatesView.vue";
import PipelinesView from "@/views/group/project/data/workspace/PipelinesView.vue";
import ResultsView from "@/views/group/project/data/workspace/ResultsView.vue";

@Component({
  components: { ResultsView, PipelinesView, GatesView, DatasetsView },
})
export default class DataWorkspaceView extends Vue {
  tab = 0;
}
</script>
