<template>
  <canvas ref="canvas" :width="canvasWidth" :height="canvasHeight" />
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import createScatterplot, { createTextureFromUrl, createRegl } from "regl-scatterplot";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { interpolateCool, scaleSequential, rgb } from "d3";
import { experimentModule } from "@/modules/experiment";
import { IChannel } from "@/modules/experiment/models";
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import {transformCoords} from "@/utils/webglUtils";

@Component
export default class ImageViewer extends Vue {
  readonly analysisContext = analysisModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  @Prop(Number) canvasWidth;
  @Prop(Number) canvasHeight;

  points: any[] = [];
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

  get centroidsData() {
    return this.analysisContext.getters.centroidsData;
  }

  // @Watch("selectedChannels")
  // onSelectedChannelsChanged(channels: IChannel[]) {
  //   if (channels && channels.length > 0) {
  //     this.experimentContext.actions.getChannelStackImage();
  //     if (this.activeDataset && this.activeAcquisitionId) {
  //       this.analysisContext.actions.getCentroidsData({
  //         dataset_id: this.activeDataset.id,
  //         acquisition_id: this.activeAcquisitionId,
  //       });
  //     }
  //   } else {
  //     this.experimentContext.mutations.setChannelStackImage(null);
  //   }
  // }

  @Watch("channelStackImage")
  onChannelStackImageChanged(value) {
    if (value) {
      this.update("channelStackImage");
    }
  }

  @Watch("centroidsData")
  onCentroidsDataChanged(value) {
    if (value) {
      this.update("centroidsData");
    }
  }

  private update(source) {
    if (this.channelStackImage) {
      const regl = this.scatterplot.get("regl");

      const img = new Image();
      img.crossOrigin = "";
      img.src = this.channelStackImage as any;
      img.onload = () => {
        this.scatterplot.set({
          backgroundImage: regl.texture(img),
        });

        if (this.centroidsData) {
          const x = transformCoords(this.centroidsData, 800, 800);
          console.log(x)

          this.points = new Array(1000).fill(0).map(() => [
            -1 + Math.random() * 2, // x
            -1 + Math.random() * 2, // y
            Math.random(), // category
            "1_111", // value
          ]);
          this.scatterplot.draw(this.points);
        } else {
         this.scatterplot.draw([]);
        }
      };
    }
  }

  pointoverHandler(pointId) {
    const [x, y, category, value] = this.points[pointId];
    // console.log(`X: ${x}\nY: ${y}\nCategory: ${category}\nValue: ${value}`);
  }

  pointoutHandler(pointId) {
    const [x, y, category, value] = this.points[pointId];
    // console.log(`X: ${x}\nY: ${y}\nCategory: ${category}\nValue: ${value}`);
  }

  selectHandler({ points: selectedPoints }) {
    console.log("Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const cell_ids: number[] = [];
      for (const i of this.selection) {
        const [acquisition_id, cell_id] = this.points[i][3].split("_");
        if (this.activeAcquisitionId === Number(acquisition_id)) {
          cell_ids.push(Number(cell_id));
        }
      }
      if (this.applyMask) {
        // this.settingsContext.mutations.setMaskSettings({
        //   ...this.settingsContext.getters.maskSettings,
        //   cell_ids: value,
        // });
        // this.experimentContext.actions.getGatedMaskImage(cell_ids);
        // this.experimentContext.actions.getChannelStackImage();
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

    this.scatterplot.subscribe("pointover", this.pointoverHandler);
    this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);

    window.addEventListener("resize", this.resizeHandler);
  }

  async mounted() {
    this.initViewer();
  }

  beforeDestroy() {
    window.removeEventListener("resize", this.resizeHandler);
    if (this.scatterplot) {
      this.scatterplot.destroy();
    }
  }
}
</script>
