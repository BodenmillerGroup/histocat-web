<template>
  <LoadingView v-if="!experimentData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <DataWorkspaceView />
    </div>
    <AnalysisView />
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import AnalysisView from "@/views/main/experiment/data/analysis/AnalysisView.vue";
import DataWorkspaceView from "@/views/main/experiment/data/workspace/DataWorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    AnalysisView,
    DataWorkspaceView,
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

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get layoutClass() {
    if (!this.showWorkspace) {
      return "layout-without-workspace py-0";
    }
    return "layout-full py-0";
  }
}
</script>

<style scoped>
.layout-full {
  display: grid;
  grid-template-columns: 380px 1fr;
  grid-template-rows: auto;
}
.layout-without-workspace {
  display: grid;
  grid-template-columns: 1fr;
  grid-template-rows: auto;
}
</style>
