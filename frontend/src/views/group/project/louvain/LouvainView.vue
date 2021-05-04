<template>
  <v-row v-if="activeResult && activeResult.output.louvain">
    <v-col cols="6">
      <PlotView title="Louvain (violin)" :src="`${url}?plot_name=stacked_violin_louvain`" />
    </v-col>
    <v-col cols="6">
      <PlotView title="Louvain (dotplot)" :src="`${url}?plot_name=dotplot_louvain`" />
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import PlotView from "@/components/PlotView.vue";
import { apiUrl } from "@/env";
import { cellsModule } from "@/modules/cells";

@Component({
  components: {
    PlotView,
  },
})
export default class LouvainView extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get activeResultId() {
    return this.cellsContext.getters.activeResultId;
  }

  get activeResult() {
    return this.cellsContext.getters.activeResult;
  }

  get url() {
    return `${apiUrl}/results/${this.activeResultId}/plot`;
  }
}
</script>
