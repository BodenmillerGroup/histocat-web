<template>
  <div id="canvasContainer">
    <canvas :id="canvas2d" :ref="canvas2d"></canvas>
    <canvas :id="canvasWebGl" :ref="canvasWebGl" v-intersect="onIntersect" v-resize="onResize" />
  </div>
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Vue, Watch } from "vue-property-decorator";
import { experimentModule } from "@/modules/experiment";
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { transformCoords } from "@/utils/webglUtils";
import { centroidsModule } from "@/modules/centroids";
import { CellPoint } from "@/data/CellPoint";
import { selectionModule } from "@/modules/selection";
import { SelectedCell } from "@/modules/selection/models";
import { mainModule } from "@/modules/main";

@Component
export default class ImageViewer extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly centroidsContext = centroidsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  private readonly canvas2d = "canvas2d";
  private readonly canvasWebGl = "canvasWebGl";

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

  @Watch("selectedChannels")
  onSelectedChannelsChanged(value) {
    this.drawLegend();
  }

  @Watch("channelStackImage")
  async onChannelStackImageChanged(value) {
    if (value && this.activeAcquisitionId) {
      const img = new Image();
      img.crossOrigin = "";
      img.onload = () => {
        const regl = this.scatterplot.get("regl");
        const prevBackgroundImage = this.scatterplot.get("backgroundImage");
        this.scatterplot.set({
          backgroundImage: regl.texture(img),
        });
        if (prevBackgroundImage) {
          prevBackgroundImage.destroy();
        }

        if (this.applyMask && this.centroids && this.centroids.has(this.activeAcquisitionId!)) {
          this.points = this.centroids.get(this.activeAcquisitionId!)!;
          const x = transformCoords(this.points, this.activeAcquisition!.max_x, this.activeAcquisition!.max_y);
          this.scatterplot.draw(x);
        } else {
          this.scatterplot.draw([]);
        }
        this.scatterplot.deselect({ preventEvent: true });
      };
      img.src = value;
    }
  }

  pointoverHandler(idx: number) {
    const point = this.points[idx];
    console.log(
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nObjectNumber: ${point.objectNumber}\nValue: ${point.value}`
    );
  }

  pointoutHandler(idx: number) {
    const point = this.points[idx];
    console.log(
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nObjectNumber: ${point.objectNumber}\nValue: ${point.value}`
    );
  }

  selectHandler({ points: selectedPoints }) {
    console.log("ImageViewer Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const newSelectedCells: SelectedCell[] = [];
      for (const i of this.selection) {
        const point = this.points[i];
        const acquisitionId = point.acquisitionId;
        newSelectedCells.push(Object.freeze(new SelectedCell(acquisitionId, point.cellId, point.objectNumber)));
      }
      this.selectionContext.actions.setSelectedCells(newSelectedCells);
      if (this.applyMask) {
        this.experimentContext.actions.getChannelStackImage();
      }
    } else {
      this.selectionContext.actions.setSelectedCells([]);
      if (this.applyMask) {
        this.experimentContext.actions.getChannelStackImage();
      }
    }
  }

  deselectHandler() {
    this.selection = [];
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
    const textHeight = this.settingsContext.getters.legend.fontScale;
    ctx.font = `${textHeight}pt sans-serif`;
    let maxTextWidth = 0;
    this.selectedChannels.forEach((v, i) => {
      const channelSettings = this.activeAcquisitionId
        ? this.settingsContext!.getters.getChannelSettings(this.activeAcquisitionId, v.name)
        : undefined;
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
      const color = this.settingsContext.getters.colorMap[v.name]
        ? this.settingsContext.getters.colorMap[v.name]
        : "#ffffff";
      const channelSettings = this.activeAcquisitionId
        ? this.settingsContext!.getters.getChannelSettings(this.activeAcquisitionId, v.name)
        : undefined;
      const text = channelSettings && channelSettings.customLabel ? channelSettings.customLabel : v.label;
      ctx.fillStyle = color!;
      ctx.fillText(text, 10, (textHeight + 10) * (i + 1) + 5);
    });
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

    this.scatterplot.subscribe("pointover", this.pointoverHandler);
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
#canvas2d {
  pointer-events: none;
  position: absolute;
  z-index: 2;
}
</style>
