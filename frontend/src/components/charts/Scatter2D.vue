<template>
  <div id="canvasContainer">
    <canvas id="canvasWebGl" ref="canvasWebGl" v-intersect="onIntersect" v-resize="onResize" />
  </div>
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import {
  scaleSequential,
  rgb,
  schemeCategory10,
  interpolateReds,
} from "d3";
import { experimentModule } from "@/modules/experiment";
import { selectionModule } from "@/modules/selection";
import { equals, uniq } from "rambda";
import { CellPoint } from "@/data/CellPoint";
import { SelectedCell } from "@/modules/selection/models";
import {mainModule} from "@/modules/main";

@Component
export default class Scatter2D extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;

  points: CellPoint[] = [];
  scatterplot: any;
  selection: any[] = [];

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

  onIntersect(entries, observer, isIntersecting) {
    if (isIntersecting) {
      const canvas = this.$refs.canvasWebGl as Element;
      const { width, height } = canvas.getBoundingClientRect();
      this.scatterplot.set({ width, height });
    }
  }

  onResize() {
    this.refresh();
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
    if (!this.scatterplot) {
      return;
    }
    const canvas = this.$refs.canvasWebGl as Element;
    const { width, height } = canvas.getBoundingClientRect();
    this.scatterplot.set({ width, height });
    this.scatterplot.refresh();
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
        this.scatterplot.select(indices, { preventEvent: true });
        this.selection = indices;
      }
    }
  }

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      this.scatterplot.deselect();

      if (data.heatmap) {
        // Use heatmap value
        const values = data.heatmap!.data;
        const max = Math.max(...values);

        const normalizedValues = values.map((v) => v / max);

        this.points = data.x.data.map(
          (x, i) => new CellPoint(data.acquisitionIds[i], data.cellIds[i], x, data.y.data[i], normalizedValues[i])
        );

        // Use continuous color scale
        const colorScale = scaleSequential(interpolateReds).domain([0, values.length]);
        const colors: string[] = [];
        for (let i = 0; i < values.length; i++) {
          colors.push(rgb(colorScale(i)).hex());
        }

        this.scatterplot.set({
          colorBy: "value",
          pointColor: colors,
        });
      } else {
        // Use acquisitionId as category
        const colorMap = new Map<number, number>();
        const uniqueAcquisitionIds = uniq(data.acquisitionIds);
        uniqueAcquisitionIds.forEach((v, i) => {
          colorMap.set(v, i);
        });
        this.points = data.x.data.map(
          (x, i) =>
            new CellPoint(
              data.acquisitionIds[i],
              data.cellIds[i],
              x,
              data.y.data[i],
              colorMap.get(data.acquisitionIds[i])!
            )
        );

        // Use discrete categorical color scale
        this.scatterplot.set({
          colorBy: "category",
          pointColor: schemeCategory10,
        });
      }

      const p = this.points.map((item) => [item.x, item.y, item.value, item.value]);
      this.scatterplot.draw(p);
    }
  }

  pointoverHandler(idx: number) {
    const point = this.points[idx];
    console.log(
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nValue: ${point.value}`
    );
  }

  pointoutHandler(idx: number) {
    const point = this.points[idx];
    console.log(
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nValue: ${point.value}`
    );
  }

  selectHandler({ points: selectedPoints }) {
    console.log("Scatter2D Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const newSelectedCells = new Map<number, SelectedCell[]>();
      for (const i of this.selection) {
        const point = this.points[i];
        const acquisitionId = point.acquisitionId;
        if (!newSelectedCells.has(acquisitionId)) {
          newSelectedCells.set(acquisitionId, []);
        }
        const ids = newSelectedCells.get(acquisitionId);
        ids!.push(Object.freeze(new SelectedCell(i, point.cellId)));
      }
      this.selectionContext.mutations.setSelectedCells(newSelectedCells);
      if (this.applyMask) {
        console.log("Scatter2D getChannelStackImage");
        this.experimentContext.actions.getChannelStackImage();
      }
    } else {
      this.selectionContext.mutations.setSelectedCells(null);
      if (this.applyMask) {
        console.log("Scatter2D getChannelStackImage");
        this.experimentContext.actions.getChannelStackImage();
      }
    }
  }

  deselectHandler() {
    console.log("Deselected:", this.selection);
    this.selection = [];
  }

  private initViewer() {
    const canvas = this.$refs.canvasWebGl as Element;
    const { width, height } = canvas.getBoundingClientRect();

    this.scatterplot = createScatterplot({
      canvas,
      width,
      height,
      opacity: 1,
      pointSize: 2,
      pointSizeSelected: 1,
      pointOutlineWidth: 1,
      lassoMinDelay: 15,
    });

    // this.scatterplot.subscribe("pointover", this.pointoverHandler);
    // this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);
  }

  mounted() {
    this.initViewer();
  }

  beforeDestroy() {
    if (this.scatterplot) {
      this.scatterplot.destroy();
    }
  }
}
</script>

<style scoped>
#canvasContainer {
  height: calc(100vh - 134px);
  position: relative;
  width: 100%;
}
#canvasWebGl {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>
