<template>
  <div>
    <canvas id="canvas2d" ref="canvas2d" :width="canvasWidth" :height="canvasHeight"></canvas>
    <canvas id="canvasWebGl" ref="canvasWebGl" :width="canvasWidth" :height="canvasHeight" v-intersect="onIntersect" />
  </div>
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
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

  get showLegend() {
    return this.settingsContext.getters.legend.apply;
  }

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

  get selectedChannels() {
    return this.experimentContext.getters.selectedChannels;
  }

  onIntersect(entries, observer, isIntersecting) {
    if (isIntersecting) {
      const canvas = this.$refs.canvasWebGl as Element;
      const { width, height } = canvas.getBoundingClientRect();
      this.scatterplot.set({ width, height });
    }
  }

  @Watch("selectedChannels")
  onSelectedChannelsChanged(value) {
    this.drawLegend();
  }

  @Watch("channelStackImage")
  async onChannelStackImageChanged(value) {
    console.log("ImageViewer update");
    if (value && this.activeAcquisitionId) {
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
    const canvas = this.$refs.canvasWebGl as any;
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

  @Watch("showLegend")
  showLegendChanged(value: boolean) {
    this.drawLegend();
  }

  private drawLegend() {
    const canvas = this.$refs.canvas2d as HTMLCanvasElement;
    const ctx = canvas.getContext("2d")!;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (!this.showLegend) {
      return;
    }
    const showIntensity = this.settingsContext.getters.legend.showIntensity;
    const textHeight = this.settingsContext.getters.legend.fontScale;
    ctx.font = `${textHeight}pt sans-serif`;
    let maxTextWidth = 0;
    this.selectedChannels.forEach((v, i) => {
      const channelSettings = this.settingsContext!.getters.getChannelSettings(this.activeAcquisitionId, v.name);
      const text = channelSettings && channelSettings.customLabel ? channelSettings.customLabel : v.label;
      const textWidth = ctx.measureText(text).width;
      if (textWidth > maxTextWidth) {
        maxTextWidth = textWidth;
      }
    });

    // Draw legend rectangle
    if (this.selectedChannels.length > 0) {
      ctx.fillStyle = "rgba(0, 0, 0, 0.8)";
      ctx.fillRect(5, 5, maxTextWidth + 10, (textHeight + 10) * this.selectedChannels.length + 10);
    }

    this.selectedChannels.forEach((v, i) => {
      const color = this.settingsContext.getters.metalColorMap.get(v.name)
        ? this.settingsContext.getters.metalColorMap.get(v.name)
        : "#ffffff";
      const channelSettings = this.settingsContext!.getters.getChannelSettings(this.activeAcquisitionId, v.name);
      const text = channelSettings && channelSettings.customLabel ? channelSettings.customLabel : v.label;
      ctx.fillStyle = color!;
      ctx.fillText(text, 10, (textHeight + 10) * (i + 1) + 5);
    });
  }

  private drawScalebar() {
    const canvas = this.$refs.canvas2d as HTMLCanvasElement;
    const { width, height } = canvas.getBoundingClientRect();
    const ctx = canvas.getContext("2d")!;
    ctx.fillStyle = "rgba(255, 255, 255, 0.8)";
    ctx.fillRect(5, height - 100, 200, 3);

    // def draw_scalebar(image: np.ndarray, scalebar: ScalebarDto):
    //   height, width, _ = image.shape
    //   length = 64
    //   cv2.line(
    //       image, (width - 60, height - 60), (width - 60 - length, height - 60), (255, 255, 255), 2, cv2.LINE_4,
    //   )
    //
    //   scale_text = length
    //   if scalebar.settings is not None and "scale" in scalebar.settings:
    //       scale = scalebar.settings.get("scale")
    //       if scale is not None and scale != "":
    //           scale_text = int(length * float(scale))
    //   cv2.putText(
    //       image,
    //       f"{scale_text} um",
    //       (width - 60 - length, height - 30),
    //       cv2.FONT_HERSHEY_PLAIN,
    //       1,
    //       (255, 255, 255),
    //       1,
    //       cv2.LINE_4,
    //   )
    //   return image
  }

  private initViewer() {
    const canvas = this.$refs.canvasWebGl as any;

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

<style scoped>
#canvasWebGl,
#canvas2d {
  position: absolute;
}
#canvas2d {
  z-index: 2;
  pointer-events: none;
}
</style>
