<template>
  <div :id="plotId" :ref="plotId" class="scatter-plot-2d" v-intersect="onIntersect">
    <!-- Plotly chart will be drawn inside this DIV -->
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import Plotly from "plotly.js-dist-min";
import { cellsModule } from "@/modules/cells";
import { ICell } from "@/modules/cells/models";
import { uiModule } from "@/modules/ui";
import { isEqual } from "lodash-es";

const pointSize = 3;
const unselectedPointSize = 2;

@Component
export default class ScatterPlot2d extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly uiContext = uiModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  @Prop({ type: String, required: true }) readonly plotId;
  @Prop({ type: String, required: true }) readonly title;
  @Prop({ type: String, required: false }) readonly xAxisTitle;
  @Prop({ type: String, required: false }) readonly yAxisTitle;
  @Prop({ type: String, required: true }) readonly mapping!: string;
  @Prop() readonly data!: Readonly<Map<number, ICell[]>>;
  @Prop({ type: Boolean, required: true }) readonly ignoreSelection!: boolean;

  get applyMask() {
    return this.uiContext.getters.showMask;
  }

  get heatmap() {
    return this.cellsContext.getters.heatmap;
  }

  get selectedCells() {
    return this.cellsContext.getters.selectedCells;
  }

  onIntersect(entries, observer, isIntersecting) {
    if (isIntersecting) {
      this.refresh();
    }
  }

  public refresh() {
    try {
      const plot = this.$refs[this.plotId] as Element;
      if (plot) {
        const { width, height } = plot.getBoundingClientRect();
        Plotly.relayout(plot, {
          width,
          height,
        });
      }
    } catch (e) {
      // TODO: find more elegant way to avoid exception during component dragging
    }
  }

  private getTraces(data: Readonly<Map<number, ICell[]>>) {
    const traces: any[] = [];
    data.forEach((cells, k) => {
      const filteredCells = cells.filter((v) => v.mappings[this.mapping]);
      traces.push({
        type: "scattergl",
        mode: "markers",
        name: `Acquisition ${k}`,
        x: filteredCells.map((v) => v.mappings[this.mapping][0]),
        y: filteredCells.map((v) => v.mappings[this.mapping][1]),
        text: filteredCells.map((v) => `CellID: ${v.cellId}`),
        customdata: filteredCells,
        marker: {
          size: pointSize,
          color: filteredCells.map((v) => v.color),
        },
        selected: {
          marker: {
            size: pointSize,
            opacity: 1.0,
            color: "#4151b6",
          },
        },
        unselected: {
          marker: {
            size: unselectedPointSize,
            opacity: 0.1,
          },
        },
      });
    });
    return traces;
  }

  private refreshOnDataChange(data: Readonly<Map<number, ICell[]>>) {
    const traces = this.getTraces(data);
    const layout = {
      showlegend: true,
      hovermode: "closest",
      dragmode: "lasso",
      autosize: true,
      xaxis: {
        title: this.xAxisTitle,
        spikesnap: "cursor",
        spikemode: "across",
        spikethickness: 1,
        spikedash: "solid",
        showspikes: true,
        spikecolor: "grey",
      },
      yaxis: {
        title: this.yAxisTitle,
        spikesnap: "cursor",
        spikemode: "across",
        spikethickness: 1,
        spikedash: "solid",
        showspikes: true,
        spikecolor: "grey",
      },
    };

    try {
      const plot = this.$refs[this.plotId] as any;
      Plotly.react(plot, traces, layout);
    } catch (e) {
      // TODO: find more elegant way to avoid exception during component dragging
    }
  }

  @Watch("data")
  dataChanged(data: Map<number, ICell[]>) {
    this.refreshOnDataChange(data);
  }

  @Watch("selectedCells")
  selectedCellsChanged(selectedCells: ICell[], oldSelectedCells: ICell[]) {
    if (!isEqual(selectedCells, oldSelectedCells)) {
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

      try {
        const plot = this.$refs[this.plotId] as any;
        Plotly.update(plot, updatedData);
      } catch (e) {
        // TODO: find more elegant way to avoid exception during component dragging
      }
    }
  }

  private initPlot() {
    const traces = this.data.size > 0 ? this.getTraces(this.data) : [];
    const initLayout = {};
    const initConfig = {
      scrollZoom: true,
      displaylogo: false,
      displayModeBar: true,
      modeBarButtonsToRemove: ["toggleSpikelines", "hoverCompareCartesian"],
      responsive: true,
    };

    try {
      const plot = this.$refs[this.plotId] as any;

      Plotly.newPlot(plot, traces, initLayout, initConfig);

      plot.on("plotly_selected", (eventData) => {
        if (eventData) {
          if (eventData.points.length > 0) {
            const selectedCells: string[] = [];
            eventData.points.forEach((point, i) => {
              const cellPoint = point.customdata as ICell;
              selectedCells.push(cellPoint.cellId);
            });
            this.cellsContext.mutations.setSelectedCellIds(selectedCells);
            if (this.applyMask) {
              this.projectsContext.actions.getChannelStackImage();
            }
          }
        }
      });

      plot.on("plotly_deselect", () => {
        this.cellsContext.mutations.setSelectedCellIds([]);
        if (this.applyMask) {
          this.projectsContext.actions.getChannelStackImage();
        }
      });
    } catch (e) {
      // TODO: find more elegant way to avoid exception during layout refresh
    }
  }

  mounted() {
    this.initPlot();
  }

  beforeDestroy() {
    try {
      Plotly.purge(this.$refs[this.plotId]);
    } catch (e) {
      // TODO: check how plotly WebGl context should be destroyed
    }
  }
}
</script>

<style scoped>
.scatter-plot-2d {
  height: 100%;
  width: 100%;
  position: relative;
}
</style>
