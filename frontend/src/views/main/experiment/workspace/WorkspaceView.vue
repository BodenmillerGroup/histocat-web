<template>
  <v-tabs v-model="tabWorkspace">
    <v-tab>Slides</v-tab>
    <v-tab>Datasets</v-tab>
    <v-tab-item eager>
      <SlidesTreeView :experiment="experiment"/>
    </v-tab-item>
    <v-tab-item>
      <DatasetsView/>
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
  import { IExperiment } from '@/modules/experiment/models';
  import DatasetsView from '@/views/main/experiment/workspace/dataset/DatasetsView.vue';
  import SlidesTreeView from '@/views/main/experiment/workspace/SlidesTreeView.vue';
  import { Component, Prop, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      DatasetsView,
      SlidesTreeView,
    },
  })
  export default class WorkspaceView extends Vue {
    @Prop(Object) experiment!: IExperiment;

    tabWorkspace = 0;
  }
</script>
