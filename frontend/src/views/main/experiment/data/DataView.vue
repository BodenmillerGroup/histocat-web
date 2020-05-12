<template>
  <LoadingView v-if="!experimentData" text="Loading..." />
  <v-container v-else fluid class="px-1 py-0">
    <v-row no-gutters>
      <v-col v-show="showWorkspace" class="pr-1" xs="3" sm="3" md="3" lg="3" xl="2">
        <WorkspaceView :experiment="experimentData" />
      </v-col>
      <v-col :cols="viewerColumns">
        <AnalysisView />
      </v-col>
    </v-row>
  </v-container>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import AnalysisView from "@/views/main/experiment/data/analysis/AnalysisView.vue";
import WorkspaceView from "@/views/main/experiment/data/workspace/WorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    AnalysisView,
    WorkspaceView,
    LoadingView,
  },
})
export default class ExperimentView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

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
    const cols = this.$vuetify.breakpoint.name === "xl" ? 10 : 9;
    return this.showWorkspace ? cols : 12;
  }
}
</script>
