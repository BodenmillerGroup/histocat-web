<template>
  <LoadingView v-if="!experimentData" text="Loading..."/>
  <v-container
    v-else
    fluid
    grid-list-md
    pa-1
  >
    <v-layout row>
      <v-flex v-show="showWorkspace" md3>
        <WorkspaceView :experiment="experimentData"/>
      </v-flex>
      <v-flex :class="viewerClass">
        <v-tabs v-model="tabExperiment">
          <v-tab>Image</v-tab>
          <v-tab>Analysis</v-tab>
          <!--          <v-tab>Workflow</v-tab>-->
          <v-tab-item>
            <ImageView/>
          </v-tab-item>
          <v-tab-item>
            <AnalysisView/>
          </v-tab-item>
          <!--          <v-tab-item>-->
          <!--            <WorkflowTab/>-->
          <!--          </v-tab-item>-->
        </v-tabs>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script lang="ts">
  import LoadingView from '@/components/LoadingView.vue';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { WebSocketManager } from '@/WebSocketManager';
  import AnalysisView from '@/views/main/experiment/analysis/AnalysisView.vue';
  import ImageView from '@/views/main/experiment/image/ImageView.vue';
  // import WorkflowTab from '@/views/main/experiment/workflow/WorkflowTab.vue';
  import WorkspaceView from '@/views/main/experiment/workspace/WorkspaceView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      AnalysisView,
      WorkspaceView,
      // WorkflowTab,
      ImageView,
      LoadingView,
    },
  })
  export default class ExperimentView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    tabExperiment = 0;

    get experiment() {
      return this.experimentContext.getters.activeExperiment;
    }

    get experimentData() {
      return this.experiment && this.experiment.slides ? this.experiment : undefined;
    }

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get viewerClass() {
      return this.showWorkspace ? 'md9' : 'md12';
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      this.experimentContext.mutations.setActiveExperimentId(experimentId);
      await this.experimentContext.actions.getExperimentData(experimentId);
      WebSocketManager.connect(experimentId);
    }

    beforeDestroy() {
      WebSocketManager.close();
      this.experimentContext.mutations.resetExperiment();
    }
  }
</script>
