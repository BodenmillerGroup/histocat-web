<template>
  <v-tabs v-model="tab">
    <v-tab>Slides</v-tab>
    <v-tab>Datasets</v-tab>
    <v-tab-item>
      <SlidesTreeView :experiment="experiment" />
    </v-tab-item>
    <v-tab-item>
      <DatasetsView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import { IExperiment } from "@/modules/experiment/models";
import DatasetsView from "@/views/main/experiment/workspace/dataset/DatasetsView.vue";
import SlidesTreeView from "@/views/main/experiment/workspace/slide/SlidesTreeView.vue";
import { Component, Prop, Vue } from "vue-property-decorator";

@Component({
  components: {
    DatasetsView,
    SlidesTreeView
  }
})
export default class WorkspaceView extends Vue {
  @Prop(Object) readonly experiment!: IExperiment;

  tab = 0;
}
</script>
