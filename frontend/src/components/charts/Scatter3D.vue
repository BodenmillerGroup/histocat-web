<template>
  <div id="scatter3dContainer">
    <div ref="scatter3d" id="scatter3d"><!-- Plotly chart will be drawn inside this DIV --></div>
  </div>
</template>

<script lang="ts">
import { IChart3DData } from "@/modules/analysis/models";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import Plotly from "plotly.js/dist/plotly";
import { mainModule } from "@/modules/main";

@Component
export default class Scatter3D extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  @Watch("showWorkspace")
  showWorkspaceChanged(value) {
    this.refresh();
  }

  @Watch("showOptions")
  showOptionsChanged(value) {
    this.refresh();
  }

  refresh() {
    const canvasElement = this.$refs.scatter3d as Element;
    const { width, height } = canvasElement.getBoundingClientRect();
    Plotly.relayout("scatter3d", {
      width,
      height,
    });
  }

  @Watch("data")
  dataChanged(data: IChart3DData) {
    if (data) {
      const plotlyData = [
        {
          type: "scatter3d",
          mode: "markers",
          marker: {
            size: 2,
          },
          x: data.x.data,
          y: data.y.data,
          z: data.z?.data,
        },
      ];

      const layout = {
        title: this.title,
        showlegend: true,
        scene: {
          xaxis: {
            title: data.x.label,
          },
          yaxis: {
            title: data.y.label,
          },
          zaxis: {
            title: data.z!.label,
          },
        },
      };

      const config = {
        scrollZoom: true,
        displaylogo: false,
        displayModeBar: true,
        modeBarButtonsToRemove: ["toggleSpikelines"],
        responsive: true,
      };

      Plotly.react("scatter3d", plotlyData, layout, config);
    }
  }
}
</script>

<style scoped>
#scatter3dContainer {
  height: calc(100vh - 84px);
  position: relative;
  width: 100%;
}
#scatter3d {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>
