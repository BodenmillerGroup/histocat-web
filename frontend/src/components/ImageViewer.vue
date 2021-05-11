<template>
  <div id="canvasContainer" ref="canvasContainer" v-intersect="onIntersect" v-resize="onResize">
    <canvas id="canvas2d" ref="canvas2d" />
    <canvas id="canvasWebGl" ref="canvasWebGl" />
    <div id="cellTooltipContainer">
      <Tooltip id="cellTooltip" v-if="cellInfo" :cell="cellInfo" />
    </div>
  </div>
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Vue, Watch } from "vue-property-decorator";
import { projectsModule } from "@/modules/projects";
import { analysisModule } from "@/modules/analysis";
import { transformToWebGl, transformFromWebGl } from "@/utils/webglUtils";
import { IRegionStatsSubmission } from "@/modules/analysis/models";
import { cellsModule } from "@/modules/cells";
import { ICell } from "@/modules/cells/models";
import { uiModule } from "@/modules/ui";
import Tooltip from "@/components/tooltip/Tooltip.vue";

@Component({
  components: { Tooltip },
})
export default class ImageViewer extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  points: ICell[] = [];
  scatterplot: any;
  selection: any[] = [];

  cellInfo: ICell | null = null;

  get applyMask() {
    return this.uiContext.getters.maskMode === "mask";
  }

  get showLegend() {
    return this.settingsContext.getters.legend.apply;
  }

  get activeAcquisitionId() {
    return this.projectsContext.getters.activeAcquisitionId;
  }

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.projectsContext.getters.channelStackImage;
  }

  get cellsByAcquisition() {
    return this.cellsContext.getters.cellsByAcquisition;
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
    return this.uiContext.getters.mouseMode;
  }

  @Watch("mouseMode")
  mouseModeChanged(value) {
    this.scatterplot.set({
      mouseMode: value,
    });
  }

  onIntersect(entries, observer, isIntersecting) {
    if (isIntersecting) {
      const canvasContainer = this.$refs.canvasContainer as Element;
      const { width, height } = canvasContainer.getBoundingClientRect();
      this.scatterplot.set({ width, height });
    }
  }

  onResize() {
    this.refresh();
  }

  @Watch("dashboardShowDrawer")
  showWorkspaceChanged(value) {
    this.refresh();
  }

  refresh() {
    if (!this.scatterplot) {
      return;
    }
    const canvasContainer = this.$refs.canvasContainer as Element;
    const { width, height } = canvasContainer.getBoundingClientRect();
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

        if (this.applyMask && this.cellsByAcquisition && this.cellsByAcquisition.has(this.activeAcquisitionId!)) {
          this.points = this.cellsByAcquisition.get(this.activeAcquisitionId!)!;
          const x = transformToWebGl(this.points, this.activeAcquisition!.max_x, this.activeAcquisition!.max_y);
          this.scatterplot.draw(x);
        } else {
          this.scatterplot.draw([]);
        }
        // this.scatterplot.deselect({ preventEvent: true });
      };
      img.src = value;
    }
  }

  pointoverHandler(idx: number) {
    const point = this.points[idx];
    this.cellInfo = point;
  }

  pointoutHandler(idx: number) {
    const point = this.points[idx];
    this.cellInfo = null;
  }

  selectHandler({ points: selectedPoints }) {
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const selectedCells: string[] = [];
      for (const i of this.selection) {
        const point = this.points[i];
        selectedCells.push(point.cellId);
      }
      this.cellsContext.mutations.setSelectedCellIds(selectedCells);
      if (this.applyMask) {
        this.projectsContext.actions.getChannelStackImage();
      }
    } else {
      this.cellsContext.mutations.setSelectedCellIds([]);
      if (this.applyMask) {
        this.cellInfo = null;
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
      this.calculateRegionStats(
        transformFromWebGl(data.coordinates, this.activeAcquisition!.max_x, this.activeAcquisition!.max_y)
      );
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

    this.scatterplot = createScatterplot({
      syncEvents: true,
      canvas: canvas,
      opacity: 1,
      pointSize: 2,
      pointSizeSelected: 0,
      pointOutlineWidth: 0,
      pointColor: [0.66, 0.66, 0.66, 1],
      pointColorActive: [0, 0.55, 1, 1],
      pointColorHover: [1, 1, 1, 1],
      lassoMinDelay: 15,
      lassoClearEvent: "lassoEnd",
      showReticle: false,
      deselectOnDblClick: true,
      deselectOnEscape: true,
      mouseMode: "panZoom",
      keyMap: { shift: "lasso", ctrl: "merge" },
    });

    this.scatterplot.subscribe("pointover", this.pointoverHandler);
    this.scatterplot.subscribe("pointout", this.pointoutHandler);
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
  height: 100%;
  width: 100%;
  position: absolute;
}
#canvasWebGl {
  height: 100%;
  width: 100%;
  position: absolute;
}
#canvas2d {
  pointer-events: none;
  position: absolute;
  z-index: 2;
}
#cellTooltipContainer {
  position: relative;
}
#cellTooltip {
  margin: 4px;
  background: black;
  color: white;
  opacity: 0.8;
  pointer-events: none;
  position: absolute;
  z-index: 2;
  top: 0;
  right: 0;
  height: auto;
  width: auto;
  font-size: 80%;
  text-align: right;
}
</style>
