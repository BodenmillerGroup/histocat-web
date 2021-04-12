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
import { projectsModule } from "@/modules/projects";
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { transformToWebGl, transformFromWebGl } from "@/utils/webglUtils";
import { centroidsModule } from "@/modules/centroids";
import { mainModule } from "@/modules/main";
import { resultsModule } from "@/modules/results";
import { ICellPoint, ISelectedCell } from "@/modules/results/models";
import { IRegionStatsSubmission } from "@/modules/analysis/models";

@Component
export default class ImageViewer extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly centroidsContext = centroidsModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);

  private readonly canvas2d = "canvas2d";
  private readonly canvasWebGl = "canvasWebGl";

  points: ICellPoint[] = [];
  scatterplot: any;
  selection: any[] = [];

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.mode === "mask";
  }

  get showLegend() {
    return this.settingsContext.getters.legend.apply;
  }

  get activeAcquisitionId() {
    return this.projectsContext.getters.activeAcquisitionId;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.projectsContext.getters.channelStackImage;
  }

  get centroids() {
    return this.centroidsContext.getters.centroids;
  }

  get selectedChannels() {
    return this.projectsContext.getters.selectedChannels;
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  get regionsEnabled() {
    return this.analysisContext.getters.regionsEnabled;
  }

  get mouseMode() {
    return this.settingsContext.getters.mouseMode;
  }

  @Watch("mouseMode")
  mouseModeChanged(value) {
    this.scatterplot.set({
      mouseMode: value,
    });
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
        const prevBackgroundImage = this.scatterplot.get("backgroundImage");
        const regl = this.scatterplot.get("regl");
        this.scatterplot.set({
          backgroundImage: regl.texture(img),
          aspectRatio: this.activeAcquisition!.max_x / this.activeAcquisition!.max_y,
        });
        if (prevBackgroundImage) {
          prevBackgroundImage.destroy();
        }

        if (this.applyMask && this.centroids && this.centroids.has(this.activeAcquisitionId!)) {
          this.points = this.centroids.get(this.activeAcquisitionId!)!;
          const x = transformToWebGl(this.points, this.activeAcquisition!.max_x, this.activeAcquisition!.max_y);
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
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nObjectNumber: ${point.objectNumber}\nColor: ${point.color}`
    );
  }

  pointoutHandler(idx: number) {
    const point = this.points[idx];
    console.log(
      `X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nObjectNumber: ${point.objectNumber}\nColor: ${point.color}`
    );
  }

  selectHandler({ points: selectedPoints }) {
    console.log("ImageViewer Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const newSelectedCells: ISelectedCell[] = [];
      for (const i of this.selection) {
        const point = this.points[i];
        const acquisitionId = point.acquisitionId;
        newSelectedCells.push(
          Object.freeze({ acquisitionId: acquisitionId, cellId: point.cellId, objectNumber: point.objectNumber })
        );
      }
      this.resultsContext.mutations.setSelectedCells(newSelectedCells);
      if (this.applyMask) {
        this.projectsContext.actions.getChannelStackImage();
      }
    } else {
      this.resultsContext.mutations.setSelectedCells([]);
      if (this.applyMask) {
        this.projectsContext.actions.getChannelStackImage();
      }
    }
  }

  deselectHandler() {
    this.selection = [];
  }

  calculateRegionStats(coordinates: number[][]) {
    if (!this.activeProjectId || !this.activeAcquisitionId) {
      return;
    }
    const params: IRegionStatsSubmission = {
      project_id: this.activeProjectId,
      acquisition_id: this.activeAcquisitionId,
      region_polygon: coordinates,
    };
    return this.analysisContext.actions.calculateRegionStats(params);
  }

  lassoEndHandler(data: { coordinates: [number, number][] }) {
    if (this.regionsEnabled) {
      this.calculateRegionStats(transformFromWebGl(data.coordinates, 800, 800));
    }
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
      const text = v.customLabel;
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
      const color = this.settingsContext.getters.channelsSettings[v.name]
        ? this.settingsContext.getters.channelsSettings[v.name].color
        : "#ffffff";
      const text = v.customLabel;
      ctx.fillStyle = color!;
      ctx.fillText(text, 10, (textHeight + 10) * (i + 1) + 5);
    });
  }

  private initViewer() {
    const canvas = this.$refs.canvasWebGl as Element;
    const { width, height } = canvas.getBoundingClientRect();

    this.scatterplot = createScatterplot({
      syncEvents: true,
      canvas: canvas,
      width: width,
      height: height,
      opacity: 1,
      pointSize: 2,
      pointSizeSelected: 1,
      pointOutlineWidth: 1,
      pointColor: [255, 140, 0],
      lassoMinDelay: 15,
      lassoClearEvent: "lassoEnd",
      showReticle: false,
      deselectOnDblClick: true,
      deselectOnEscape: true,
      mouseMode: "panZoom",
    });

    this.scatterplot.subscribe("pointover", this.pointoverHandler);
    // this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);
    this.scatterplot.subscribe("lassoEnd", this.lassoEndHandler);
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
