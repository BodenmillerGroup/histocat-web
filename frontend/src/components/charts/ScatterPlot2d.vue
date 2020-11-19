<template>
  <div :id="plotId" :ref="plotId" class="plot" v-intersect="onIntersect" v-resize="onResize">
    <!-- Plotly chart will be drawn inside this DIV -->
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { settingsModule } from "@/modules/settings";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import Plotly from "plotly.js/dist/plotly";
import { ICellPoint, ISelectedCell } from "@/modules/results/models";
import { resultsModule } from "@/modules/results";

@Component
export default class ScatterPlot2d extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);

  @Prop({ type: String, required: true }) plotId;
  @Prop({ type: String, required: true }) title;
  @Prop({ type: String, required: true }) xAxisTitle;
  @Prop({ type: String, required: true }) yAxisTitle;
  @Prop({ type: Map, required: true }) data!: Map<number, ICellPoint[]>;
  @Prop({ type: Boolean, required: true }) ignoreSelection!: boolean;

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  get heatmap() {
    return this.resultsContext.getters.heatmap;
  }

  get selectedCells() {
    return this.resultsContext.getters.selectedCells;
  }

  onIntersect(entries, observer, isIntersecting) {
    if (isIntersecting) {
      this.refresh();
    }
  }

  onResize() {
    // this.refresh();
  }

  refresh() {
    const plotElement = this.$refs[this.plotId] as Element;
    const { width, height } = plotElement.getBoundingClientRect();
    Plotly.relayout(this.plotId, {
      width,
      height,
    });
  }

  private refreshOnDataChange(data: Map<number, ICellPoint[]>) {
    if (data) {
      const traces: any[] = [];
      data.forEach((v, k) => {
        traces.push({
          type: "scattergl",
          mode: "markers",
          name: `Acquisition ${k}`,
          x: v.map((v) => v.x),
          y: v.map((v) => v.y),
          text: v.map((v) => `CellID: ${v.cellId}`),
          customdata: v,
          marker: this.heatmap
            ? {
                size: 3,
                color: v.map((v) => v.color),
                colorscale: "Jet",
              }
            : {
                size: 3,
              },
          unselected: {
            marker: {
              opacity: 0.1,
            },
          },
        });
      });

      const layout = {
        showlegend: true,
        xaxis: {
          title: this.xAxisTitle,
        },
        yaxis: {
          title: this.yAxisTitle,
        },
        hovermode: "closest",
        dragmode: "lasso",
        autosize: true,
      };

      Plotly.react(this.plotId, traces, layout);
    } else {
      Plotly.react(this.plotId, [], {});
    }
  }

  @Watch("data")
  dataChanged(data: Map<number, ICellPoint[]>) {
    this.refreshOnDataChange(data);
  }

  @Watch("selectedCells")
  selectedCellsChanged(selectedCells: ISelectedCell[]) {
    if (this.data) {
      const selectedpoints: any[] = [];
      this.data.forEach((v, k) => {
        // TODO: Important!! ObjectNumber starts from 1, so index should be ObjectNumber - 1
        selectedpoints.push(
          selectedCells.length > 0
            ? selectedCells.filter((v) => v.acquisitionId === k).map((v) => v.objectNumber - 1)
            : null
        );
      });

      const updatedData: any = {
        selectedpoints: selectedpoints,
      };

      if (this.ignoreSelection) {
        updatedData.unselected = {
          marker: {
            opacity: 0,
          },
        };
      }

      Plotly.update(this.plotId, updatedData);
    }
  }

  private initPlot() {
    const initData = [];
    const initLayout = {};
    const initConfig = {
      scrollZoom: true,
      displaylogo: false,
      displayModeBar: true,
      modeBarButtonsToRemove: ["toggleSpikelines", "hoverCompareCartesian"],
      responsive: true,
    };

    Plotly.react(this.plotId, initData, initLayout, initConfig);

    const plot = this.$refs[this.plotId] as any;
    plot.on("plotly_selected", (eventData) => {
      if (eventData) {
        if (eventData.points.length > 0) {
          // console.log(eventData.points);
          const newSelectedCells: ISelectedCell[] = [];
          eventData.points.forEach((point, i) => {
            const cellPoint = point.customdata as ICellPoint;
            newSelectedCells.push({
              acquisitionId: cellPoint.acquisitionId,
              cellId: cellPoint.cellId,
              objectNumber: cellPoint.objectNumber,
            });
          });
          this.resultsContext.actions.setSelectedCells(newSelectedCells);
          if (this.applyMask) {
            this.projectsContext.actions.getChannelStackImage();
          }
        }
      }
    });

    plot.on("plotly_deselect", () => {
      this.resultsContext.actions.setSelectedCells([]);
      if (this.applyMask) {
        this.projectsContext.actions.getChannelStackImage();
      }
    });

    if (this.data) {
      this.refreshOnDataChange(this.data);
    }
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
