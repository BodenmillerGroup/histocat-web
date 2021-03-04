<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <ImageWorkspaceView :projectData="projectData" />
    </div>
    <VisualizationView />
    <div v-show="showOptions">
      <OptionsView />
    </div>
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import VisualizationView from "@/views/group/project/image/visualization/VisualizationView.vue";
import ImageWorkspaceView from "@/views/group/project/image/workspace/ImageWorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";
import OptionsView from "@/views/group/project/image/options/OptionsView.vue";

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
  readonly projectsContext = projectsModule.context(this.$store);

  get projectData() {
    return this.projectsContext.getters.projectData;
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
