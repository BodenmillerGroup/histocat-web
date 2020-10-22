<template>
  <div :id="plotId" :ref="plotId" class="plot"><!-- Plotly chart will be drawn inside this DIV --></div>
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { SelectedCell } from "@/modules/selection/models";
import { selectionModule } from "@/modules/selection";
import { CellPoint } from "@/data/CellPoint";
import Plotly from "plotly.js/dist/plotly";
import { linearRegression, linearRegressionLine } from "simple-statistics";

@Component
export default class ScatterPlot2d extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(String) plotId;
  @Prop(Object) data;
  @Prop(String) title;
  @Prop(Boolean) showRegression;

  selection: SelectedCell[] = [];
  hasHeatmap = false;
  xAxisTitle = "";
  yAxisTitle = "";
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

  get selectedCells() {
    return this.selectionContext.getters.selectedCells;
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
    const plotElement = this.$refs[this.plotId] as Element;
    const { width, height } = plotElement.getBoundingClientRect();
    Plotly.relayout(this.plotId, {
      width,
      height,
    });
  }

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      this.hasHeatmap = !!data.heatmap;
      this.xAxisTitle = data.x.label;
      this.yAxisTitle = data.y.label;

      this.points = new Map<number, CellPoint[]>();
      const regressionData: number[][] = [];
      for (let i = 0; i < data.cellIds.length; i++) {
        const cellPoint = Object.freeze(
          new CellPoint(
            data.acquisitionIds[i],
            data.cellIds[i],
            data.objectNumbers[i],
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

        if (this.showRegression) {
          regressionData.push([cellPoint.x, cellPoint.y]);
        }
      }

      const traces: any[] = [];
      this.points.forEach((v, k) => {
        traces.push({
          type: "scattergl",
          mode: "markers",
          name: `Acquisition ${k}`,
          x: v.map((v) => v.x),
          y: v.map((v) => v.y),
          text: v.map((v) => `Cell ID: ${v.cellId}`),
          customdata: v,
          marker: this.hasHeatmap
            ? {
                size: 3,
                color: v.map((v) => v.value),
                colorscale: "RdBu",
              }
            : {
                size: 3,
              },
        });
      });

      if (this.showRegression) {
        const regression = linearRegression(regressionData);
        const regressionLine = linearRegressionLine(regression);
        const xMin = Math.min(...data.x.data);
        const xMax = Math.max(...data.x.data);
        traces.push({
          type: "scattergl",
          mode: "line",
          name: `Regression`,
          x: [xMin, xMax],
          y: [regressionLine(xMin), regressionLine(xMax)],
          marker: {
            size: 1,
          },
          line: {
            width: 1,
          },
        });
      }

      const layout = {
        title: this.title,
        showlegend: true,
        xaxis: {
          title: this.xAxisTitle,
        },
        yaxis: {
          title: this.yAxisTitle,
        },
        hovermode: "closest",
        autosize: true,
      };

      Plotly.react(this.plotId, traces, layout);
    }
  }

  @Watch("selectedCells")
  selectedCellsChanged(data: SelectedCell[]) {
    if (!equals(data, this.selection)) {
      this.selection = data;
      const traces: any[] = [];
      this.points.forEach((v, k) => {
        // TODO: Important!! ObjectNumber starts from 1, so index should be ObjectNumber - 1
        traces.push({
          type: "scattergl",
          mode: "markers",
          name: `Acquisition ${k}`,
          x: v.map((v) => v.x),
          y: v.map((v) => v.y),
          text: v.map((v) => `Cell ID: ${v.cellId}`),
          customdata: v,
          selectedpoints: data.filter((v) => v.acquisitionId === k).map((v) => v.objectNumber - 1),
          marker: this.hasHeatmap
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
          title: this.xAxisTitle,
        },
        yaxis: {
          title: this.yAxisTitle,
        },
      };

      Plotly.react(this.plotId, traces, layout);
    }
  }

  private initPlot() {
    const initData = [];

    const initLayout = {
      title: this.title,
      showlegend: true,
      hovermode: "closest",
    };

    const initConfig = {
      scrollZoom: true,
      displaylogo: false,
      displayModeBar: true,
      modeBarButtonsToRemove: ["toggleSpikelines"],
      responsive: true,
    };

    Plotly.react(this.plotId, initData, initLayout, initConfig);

    const plot = this.$refs[this.plotId] as any;
    plot.on("plotly_selected", (eventData) => {
      if (eventData) {
        if (eventData.points.length > 0) {
          // console.log(eventData.points);
          const newSelectedCells: SelectedCell[] = [];
          eventData.points.forEach((point, i) => {
            const cellPoint = point.customdata as CellPoint;
            newSelectedCells.push(
              Object.freeze(new SelectedCell(cellPoint.acquisitionId, cellPoint.cellId, cellPoint.objectNumber))
            );
          });
          this.selectionContext.actions.setSelectedCells(newSelectedCells);
          if (this.applyMask) {
            this.projectsContext.actions.getChannelStackImage();
          }
        } else {
          if (this.selectionContext.getters.selectedCells !== null) {
            this.selectionContext.actions.setSelectedCells([]);
            if (this.applyMask) {
              this.projectsContext.actions.getChannelStackImage();
            }
          }
        }
      }
    });
  }

  mounted() {
    this.initPlot();
  }
}
</script>

<style scoped>
.plot {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>
