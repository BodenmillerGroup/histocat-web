<template>
  <div class="widget-container">
    <ScatterPlot2d
      plot-id="pcaPlot"
      ref="pcaPlot"
      :ignore-selection="false"
      :data="cellsList"
      mapping="pca"
      title="PCA"
      x-axis-title="PC1"
      y-axis-title="PC2"
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
export default class PcaView extends Vue {
  readonly cellsContext = cellsModule.context(this.$store);

  get cellsList() {
    return this.cellsContext.getters.cellsList;
  }

  refresh() {
    if (this.$refs.pcaPlot) {
      (this.$refs.pcaPlot as any).refresh();
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
