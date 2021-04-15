<template>
  <v-tabs v-model="tab">
    <v-tab>Datasets</v-tab>
    <v-tab>Pipelines</v-tab>
    <v-tab>Gates</v-tab>
    <v-tab-item>
      <DatasetsView />
      <ResultsView />
    </v-tab-item>
    <v-tab-item>
      <PipelinesView />
    </v-tab-item>
    <v-tab-item>
      <GatesView />
      <AnnotationsView />
      <CellClassesView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import DatasetsView from "@/views/group/project/data/workspace/DatasetsView.vue";
import GatesView from "@/views/group/project/data/workspace/GatesView.vue";
import PipelinesView from "@/views/group/project/data/workspace/PipelinesView.vue";
import ResultsView from "@/views/group/project/data/workspace/ResultsView.vue";
import CellClassesView from "@/views/group/project/data/workspace/CellClassesView.vue";
import AnnotationsView from "@/views/group/project/data/workspace/AnnotationsView.vue";

@Component({
  components: { AnnotationsView, CellClassesView, ResultsView, PipelinesView, GatesView, DatasetsView },
})
export default class DataWorkspaceView extends Vue {
  tab = 0;
}
</script>
