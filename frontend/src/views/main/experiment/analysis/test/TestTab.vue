<template>
  <v-row no-gutters class="chart-container">
    <v-col :cols="columns">
      <Graph :key="graphRenderCounter" />
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Component, Vue, Watch } from "vue-property-decorator";
import Graph from "@/cellxgene/components/graph/graph.vue";
import {controlsModule} from "@/modules/controls";
import {graphModule} from "@/modules/graph";

@Component({
  components: { Graph },
})
export default class TestTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly controlsContext = controlsModule.context(this.$store);
  readonly graphContext = graphModule.context(this.$store);

  get graphRenderCounter() {
    return this.controlsContext.getters.graphRenderCounter;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    return this.showOptions ? 9 : 12;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get schema() {
    return this.graphContext.getters.schema;
  }

  get channels() {
    return this.datasetContext.getters.channels;
  }

  submit() {
    this.graphContext.actions.getSchema();
    this.graphContext.actions.getVarAnnotation("name_0");
  }

  mounted() {
    this.submit();
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
