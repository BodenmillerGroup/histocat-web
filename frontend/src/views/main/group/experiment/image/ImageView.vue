<template>
  <LoadingView v-if="!experimentData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <ImageWorkspaceView :experiment="experimentData" />
    </div>
    <VisualizationView />
    <div v-show="showOptions">
      <OptionsView />
    </div>
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import VisualizationView from "@/views/main/group/experiment/image/visualization/VisualizationView.vue";
import ImageWorkspaceView from "@/views/main/group/experiment/image/workspace/ImageWorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";
import OptionsView from "@/views/main/group/experiment/image/options/OptionsView.vue";

@Component({
  components: {
    OptionsView,
    ImageWorkspaceView,
    VisualizationView,
    LoadingView,
  },
})
export default class ImageView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

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
    if (!this.showWorkspace && this.showOptions) {
      return "layout-without-workspace py-0";
    } else if (this.showWorkspace && !this.showOptions) {
      return "layout-without-options py-0";
    } else if (!this.showWorkspace && !this.showOptions) {
      return "layout-empty py-0";
    }
    return "layout-full py-0";
  }
}
</script>

<style scoped>
.layout-full {
  display: grid;
  grid-template-columns: 380px 1fr 380px;
  grid-template-rows: auto;
}
.layout-without-workspace {
  display: grid;
  grid-template-columns: 1fr 380px;
  grid-template-rows: auto;
}
.layout-without-options {
  display: grid;
  grid-template-columns: 380px 1fr;
  grid-template-rows: auto;
}
.layout-empty {
  display: grid;
  grid-template-columns: 1fr;
  grid-template-rows: auto;
}
</style>
