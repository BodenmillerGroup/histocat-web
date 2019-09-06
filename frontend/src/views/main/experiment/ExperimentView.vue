<template>
  <LoadingView v-if="!experimentData" text="Loading..."/>
  <v-container
    v-else
    fluid
    class="px-1 py-0"
  >
    <v-row no-gutters>
      <v-col v-show="showWorkspace" cols="3" class="pr-1">
        <WorkspaceView :experiment="experimentData"/>
      </v-col>
      <v-col :cols="viewerColumns">
        <v-tabs v-model="tabExperiment">
          <v-tab>Image</v-tab>
          <v-tab>Analysis</v-tab>
          <v-tab-item>
            <ImageView/>
          </v-tab-item>
          <v-tab-item>
            <AnalysisView/>
          </v-tab-item>
        </v-tabs>
      </v-col>
    </v-row>
  </v-container>
</template>

<script lang="ts">
  import LoadingView from '@/components/LoadingView.vue';
  import { analysisModule } from '@/modules/analysis';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { WebSocketManager } from '@/utils/WebSocketManager';
  import AnalysisView from '@/views/main/experiment/analysis/AnalysisView.vue';
  import ImageView from '@/views/main/experiment/image/ImageView.vue';
  import WorkspaceView from '@/views/main/experiment/workspace/WorkspaceView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      AnalysisView,
      WorkspaceView,
      ImageView,
      LoadingView,
    },
  })
  export default class ExperimentView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);

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

    get viewerColumns() {
      return this.showWorkspace ? 9 : 12;
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      this.experimentContext.mutations.setActiveExperimentId(experimentId);
      await this.experimentContext.actions.getExperimentData(experimentId);
      WebSocketManager.connect(experimentId);
    }

    beforeDestroy() {
      WebSocketManager.close();
      if (process.env.VUE_APP_ENV !== 'development') {
        this.experimentContext.mutations.reset();
        this.analysisContext.mutations.reset();
      }
    }
  }
</script>
