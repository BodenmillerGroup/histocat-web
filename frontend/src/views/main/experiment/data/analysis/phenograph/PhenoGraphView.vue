<template>
  <div class="plotContainer">
    <div :id="plotId" ref="plot" class="plot"><!-- Plotly chart will be drawn inside this DIV --></div>
  </div>
</template>

<script lang="ts">
import { IPhenoGraphData } from "@/modules/analysis/models";
import { mainModule } from "@/modules/main";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import Plotly from "plotly.js/dist/plotly";

@Component
export default class PhenoGraphView extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  @Prop(String) plotId;
  @Prop(Object) data;
  @Prop(String) title;

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  @Watch("showWorkspace")
  showWorkspaceChanged(value: boolean) {
    this.refresh();
  }

  @Watch("showOptions")
  showOptionsChanged(value: boolean) {
    this.refresh();
  }

  refresh() {
    if (!this.data) {
      return;
    }
    const plotElement = this.$refs.plot as Element;
    const { width, height } = plotElement.getBoundingClientRect();
    Plotly.relayout(this.plotId, {
      width,
      height,
    });
  }

  @Watch("data")
  dataChanged(data: IPhenoGraphData) {
    if (data) {
      const xValues: any = data.community;
      const yValues = Object.keys(data).filter((item) => item !== "community");
      const zValues: any[] = [];
      yValues.forEach((marker) => {
        zValues.push(data[marker]);
      });

      const traces = [
        {
          type: "heatmap",
          x: xValues.map((v) => `c${v}`),
          y: yValues,
          z: zValues,
        },
      ];

      const layout = {
        title: this.title,
        showlegend: true,
        xaxis: {
          title: "Community",
        },
        yaxis: {
          title: "Marker",
        },
      };

      const config = {
        scrollZoom: true,
        displaylogo: false,
        displayModeBar: true,
        responsive: true,
        modeBarButtonsToRemove: ["toggleSpikelines"],
      };

      Plotly.react(this.plotId, traces, layout, config);
    }
  }
}
</script>

<style scoped>
.plotContainer {
  height: calc(100vh - 84px);
  position: relative;
  width: 100%;
}
.plot {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>
