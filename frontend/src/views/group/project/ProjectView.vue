<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <v-container v-else fluid class="ma-0 pa-0">
    <DataView v-show="showData" />
    <ImageView v-show="!showData" />
  </v-container>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { analysisModule } from "@/modules/analysis";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import ImageView from "@/views/group/project/image/ImageView.vue";
import DataView from "@/views/group/project/data/DataView.vue";

@Component({
  components: {
    DataView,
    ImageView,
    LoadingView,
  },
})
export default class ProjectView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  get projectData() {
    return this.projectsContext.getters.projectData;
  }

  get showData() {
    return this.mainContext.getters.showData;
  }

  async mounted() {
    const projectId = +this.$router.currentRoute.params.projectId;
    this.projectsContext.mutations.setActiveProjectId(projectId);
    await this.projectsContext.actions.getProjectData(projectId);
    WebSocketManager.connect(projectId);
  }

  beforeDestroy() {
    WebSocketManager.close();
    // if (process.env.VUE_APP_ENV !== "development") {
    //   this.$store.dispatch("reset");
    // }
    this.$store.dispatch("resetProject");
  }
}
</script>