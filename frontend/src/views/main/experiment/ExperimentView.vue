<template>
  <LoadingView v-if="!experimentData" text="Loading..." />
  <v-container v-else fluid class="px-1 py-0">
    <DataView v-show="showData" />
    <ImageView v-show="!showData" />
  </v-container>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import ImageView from "@/views/main/experiment/image/ImageView.vue";
import DataView from "@/views/main/experiment/data/DataView.vue";

@Component({
  components: {
    DataView,
    ImageView,
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

  get showData() {
    return this.mainContext.getters.showData;
  }

  async mounted() {
    const experimentId = parseInt(this.$router.currentRoute.params.experimentId, 10);
    this.experimentContext.mutations.setActiveExperimentId(experimentId);
    await this.experimentContext.actions.getExperimentData(experimentId);
    WebSocketManager.connect(experimentId);
  }

  beforeDestroy() {
    WebSocketManager.close();
    // if (process.env.VUE_APP_ENV !== "development") {
    //   this.experimentContext.mutations.reset();
    //   this.analysisContext.mutations.reset();
    // }
    this.$store.dispatch("reset");
  }
}
</script>
