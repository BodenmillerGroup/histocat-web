<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <div v-else :class="layoutClass">
    <div v-show="showWorkspace" class="pr-1">
      <DataWorkspaceView />
    </div>
    <AnalysisView />
  </div>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import AnalysisView from "@/views/group/project/data/analysis/AnalysisView.vue";
import DataWorkspaceView from "@/views/group/project/data/workspace/DataWorkspaceView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    AnalysisView,
    DataWorkspaceView,
    LoadingView,
  },
})
export default class DataView extends Vue {
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
