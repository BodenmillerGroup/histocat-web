<template>
  <v-row v-if="activeResult && activeResult.output.leiden">
    <v-col cols="6">
      <PlotView title="Leiden (violin)" :src="`${url}?plot_name=stacked_violin_leiden`" />
    </v-col>
    <v-col cols="6">
      <PlotView title="Leiden (dotplot)" :src="`${url}?plot_name=dotplot_leiden`" />
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
export default class LeidenView extends Vue {
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
