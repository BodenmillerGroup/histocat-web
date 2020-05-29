<template>
  <LoadingView v-if="!experimentData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <ImageWorkspaceView :experiment="experimentData" />
    </div>
    <div>
      <VisualizationView />
    </div>
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import VisualizationView from "@/views/main/experiment/image/visualization/VisualizationView.vue";
import ImageWorkspaceView from "@/views/main/experiment/image/workspace/ImageWorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    ImageWorkspaceView,
    VisualizationView,
    LoadingView,
  },
})
export default class ImageView extends Vue {
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

  get layoutClass() {
    return this.showWorkspace ? "layout-workspace px-1 py-0" : "px-1 py-0";
  }
}
</script>

<style scoped>
.layout-workspace {
  display: grid;
  grid-template-columns: 380px 1fr;
  grid-template-rows: auto;
}
</style>
