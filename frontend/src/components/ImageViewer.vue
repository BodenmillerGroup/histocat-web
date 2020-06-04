<template>
  <canvas ref="canvas" :width="canvasWidth" :height="canvasHeight" />
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import createScatterplot, { createTextureFromUrl, createRegl } from "regl-scatterplot/src";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { experimentModule } from "@/modules/experiment";
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { transformCoords } from "@/utils/webglUtils";
import { centroidsModule } from "@/modules/centroids";
import { CellPoint } from "@/data/CellPoint";
import { selectionModule } from "@/modules/selection";
import { SelectedCell } from "@/modules/selection/models";

@Component
export default class ImageViewer extends Vue {
  readonly analysisContext = analysisModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly centroidsContext = centroidsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Number) canvasWidth;
  @Prop(Number) canvasHeight;

  points: CellPoint[] = [];
  scatterplot: any;
  selection: any[] = [];

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  // get selectedChannels() {
  //   return this.experimentContext.getters.selectedChannels;
  // }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.experimentContext.getters.channelStackImage;
  }

  get centroids() {
    return this.centroidsContext.getters.centroids;
  }

  @Watch("channelStackImage")
  async onChannelStackImageChanged(value) {
    console.log("ImageViewer update");
    if (value && this.activeAcquisitionId) {
      const canvas = this.$refs.canvas as any;
      const { width, height } = canvas.getBoundingClientRect();
      console.log(width, height)
      // const backgroundImage = await this.scatterplot.createTextureFromUrl(
      //   `https://picsum.photos/${Math.min(640, width)}/${Math.min(640, height)}/?random`
      // );
      // this.scatterplot.set({ backgroundImage });

      const img = new Image();
      img.crossOrigin = "";
      img.onload = () => {
        const regl = this.scatterplot.get("regl");
        this.scatterplot.set({
          backgroundImage: regl.texture(img),
        });
      };
      img.src = this.channelStackImage as any;

      if (this.centroids && this.centroids.has(this.activeAcquisitionId!)) {
        this.points = this.centroids.get(this.activeAcquisitionId!)!;
        const x = transformCoords(this.points, this.activeAcquisition!.max_x, this.activeAcquisition!.max_y);
        this.scatterplot.draw(x);
      } else {
        this.scatterplot.draw([]);
      }
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
    console.log("ImageViewer Selected:", selectedPoints);
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
        console.log("ImageViewer getChannelStackImage");
        this.experimentContext.actions.getChannelStackImage();
      }
    } else {
      this.selectionContext.mutations.setSelectedCells(null);
      if (this.applyMask) {
        console.log("ImageViewer getChannelStackImage");
        this.experimentContext.actions.getChannelStackImage();
      }
    }
  }

  deselectHandler() {
    console.log("Deselected:", this.selection);
    this.selection = [];
  }

  resizeHandler() {
    const canvas = this.$refs.canvas as any;
    if (canvas) {
      const rect = canvas.getBoundingClientRect();
      if (rect) {
        this.scatterplot.set({ width: this.canvasWidth, height: this.canvasHeight });
      }
    }
  }

  refresh() {
    this.scatterplot.set({ width: this.canvasWidth, height: this.canvasHeight });
    this.scatterplot.refresh();
  }

  @Watch("canvasWidth")
  onCanvasWidthUpdate(value: number) {
    this.refresh();
  }

  private initViewer() {
    const canvas = this.$refs.canvas as any;

    this.scatterplot = createScatterplot({
      canvas,
      width: this.canvasWidth,
      height: this.canvasHeight,
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

    // window.addEventListener("resize", this.resizeHandler);
  }

  async mounted() {
    this.initViewer();
  }

  beforeDestroy() {
    // window.removeEventListener("resize", this.resizeHandler);
    if (this.scatterplot) {
      this.scatterplot.destroy();
    }
  }
}
</script>
