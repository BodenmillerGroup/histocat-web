<template>
  <div id="canvasContainer">
    <div ref="canvas" id="canvas"><!-- Plotly chart will be drawn inside this DIV --></div>
  </div>
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { interpolateReds, scaleSequential, schemeCategory10, scaleOrdinal } from "d3";
import { SelectedCell } from "@/modules/selection/models";
import { selectionModule } from "@/modules/selection";
import { CellPoint } from "@/data/CellPoint";
import Plotly from "plotly.js/dist/plotly";
import * as d3 from "d3";

@Component
export default class ScatterPlotView2 extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;

  @Prop(Boolean) showRegression;
  @Prop(String) regressionType;
  @Prop(Number) polynomialOrder;

  selection: any[] = [];
  colorScale: d3.ScaleSequential<any> | d3.ScaleOrdinal<any, any> | null = null;
  points = new Map<number, CellPoint[]>();

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get selectedCells() {
    return this.selectionContext.getters.selectedCells;
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
    const canvasElement = this.$refs.canvas as Element;
    const { width, height } = canvasElement.getBoundingClientRect();
    Plotly.relayout("canvas", {
      width,
      height,
    });
  }

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      this.points = new Map<number, CellPoint[]>();
      for (let i = 0; i < data.cellIds.length; i++) {
        const cellPoint = Object.freeze(
          new CellPoint(
            i,
            data.acquisitionIds[i],
            data.cellIds[i],
            data.x.data[i],
            data.y.data[i],
            data.heatmap ? data.heatmap.data[i] : data.acquisitionIds[i]
          )
        );
        if (!this.points.has(cellPoint.acquisitionId)) {
          this.points.set(cellPoint.acquisitionId, []);
        }
        const acquisitionPoints = this.points.get(cellPoint.acquisitionId)!;
        acquisitionPoints.push(cellPoint);
      }

      let plotlyData: any[] = [];
      this.points.forEach((v, k) => {
        plotlyData.push({
          type: "scattergl",
          mode: "markers",
          name: `Acquisition ${k}`,
          x: v.map((v) => v.x),
          y: v.map((v) => v.y),
          customdata: v,
          marker: data.heatmap
            ? {
                size: 3,
                color: v.map((v) => v.value),
              }
            : {
                size: 3,
              },
        });
      });

      const layout = {
        title: this.title,
        showlegend: true,
        xaxis: {
          title: data.x.label,
        },
        yaxis: {
          title: data.y.label,
        },
      };

      Plotly.react("canvas", plotlyData, layout);
    }
  }

  mounted() {
    const plotlyData = [
      {
        type: "scattergl",
        mode: "markers",
      },
    ];

    const layout = {
      title: this.title,
      showlegend: true,
    };

    const config = {
      scrollZoom: true,
      displaylogo: false,
      displayModeBar: true,
      modeBarButtonsToRemove: ["toggleSpikelines"],
      responsive: true,
    };

    Plotly.react("canvas", plotlyData, layout, config);

    const canvas = this.$refs.canvas as any;
    canvas.on("plotly_selected", (eventData) => {
      if (eventData) {
        if (eventData.points.length > 0) {
          const newSelectedCells = new Map<number, SelectedCell[]>();
          eventData.points.forEach((point, i) => {
            const cellPoint = point.customdata;
            if (!newSelectedCells.has(cellPoint.acquisitionId)) {
              newSelectedCells.set(cellPoint.acquisitionId, []);
            }
            const array = newSelectedCells.get(cellPoint.acquisitionId)!;
            array.push(Object.freeze(new SelectedCell(cellPoint.acquisitionId, cellPoint.index, cellPoint.cellId)));
          });
          this.selectionContext.mutations.setSelectedCells(newSelectedCells);
          if (this.applyMask) {
            this.experimentContext.actions.getChannelStackImage();
          }
        } else {
          if (this.selectionContext.getters.selectedCells !== null) {
            this.selectionContext.mutations.setSelectedCells(null);
            if (this.applyMask) {
              this.experimentContext.actions.getChannelStackImage();
            }
          }
        }
      }
    });
  }

  @Watch("selectedCells")
  selectedCellsChanged(data: Map<number, SelectedCell[]> | null) {
    if (data !== null) {
      let allCells: SelectedCell[] = [];
      data.forEach((val, key) => {
        allCells = allCells.concat(val);
      });
      const indices = allCells.map((i) => i.index);
      if (!equals(indices, this.selection)) {
        this.selection = indices;
      }
    }
  }
}
</script>

<style scoped>
#canvasContainer {
  height: calc(100vh - 84px);
  position: relative;
  width: 100%;
}
#canvas {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>
