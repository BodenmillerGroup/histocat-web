<template>
  <div class="widget-container">
    <ScatterPlot2d
      plot-id="tsnePlot"
      ref="tsnePlot"
      :ignore-selection="false"
      :data="cellsList"
      mapping="tsne"
      title="tSNE"
      x-axis-title="tSNE1"
      y-axis-title="tSNE2"
    />
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import ScatterPlot2d from "@/components/charts/ScatterPlot2d.vue";
import { cellsModule } from "@/modules/cells";

@Component({
  components: { ScatterPlot2d },
})
export default class TsneView extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get cellsList() {
    return this.cellsContext.getters.cellsList;
  }

  refresh() {
    if (this.$refs.tsnePlot) {
      (this.$refs.tsnePlot as any).refresh();
    }
  }
}
</script>

<style scoped>
.widget-container {
  width: 100%;
  height: 100%;
}
</style>
