<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <AcquisitionsView />
    </div>

    <div v-show="showOptions">
      <PanelView />
    </div>
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";
import AcquisitionsView from "@/views/group/project/segmentation/AcquisitionsView.vue";
import PanelView from "@/views/group/project/segmentation/PanelView.vue";

@Component({
  components: {
    PanelView,
    AcquisitionsView,
    LoadingView,
  },
})
export default class SegmentationView extends Vue {
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
  grid-template-columns: 400px 1fr 380px;
  grid-template-rows: auto;
}
.layout-without-workspace {
  display: grid;
  grid-template-columns: 1fr 380px;
  grid-template-rows: auto;
}
.layout-without-options {
  display: grid;
  grid-template-columns: 400px 1fr;
  grid-template-rows: auto;
}
.layout-empty {
  display: grid;
  grid-template-columns: 1fr;
  grid-template-rows: auto;
}
</style>
