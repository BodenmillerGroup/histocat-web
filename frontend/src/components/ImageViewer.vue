<template>
  <canvas ref="canvas" :width="width" :height="height" />
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { settingsModule } from "@/modules/settings";
import createScatterplot, { createTextureFromUrl, createRegl } from "regl-scatterplot";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { interpolateCool, scaleSequential, rgb } from "d3";
import { experimentModule } from "@/modules/experiment";
import { IChannel } from "@/modules/experiment/models";

@Component
export default class ImageViewer extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  @Prop(Number) width;
  @Prop(Number) height;

  points: any[] = [];
  scatterplot: any;
  selection: any[] = [];

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  get selectedChannels() {
    return this.experimentContext.getters.selectedChannels;
  }

  get metalColorMap() {
    return this.settingsContext.getters.metalColorMap;
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.experimentContext.getters.channelStackImage;
  }

  @Watch("selectedChannels")
  onSelectedChannelsChanged(channels: IChannel[]) {
    if (channels && channels.length > 0) {
      this.experimentContext.actions.getChannelStackImage();
    } else {
      this.experimentContext.mutations.setChannelStackImage(null);
    }
  }

  @Watch("channelStackImage")
  async onChannelStackImageChanged(image) {
    if (image) {
      if (this.scatterplot != null) {
        window.removeEventListener("resize", this.resizeHandler);
        this.scatterplot.destroy();
        this.scatterplot = null;
      }

      const canvas = this.$refs.canvas as any;

      const { width, height } = canvas.getBoundingClientRect();

      const regl = createRegl(canvas);
      // const backgroundImage = await createTextureFromUrl(
      //   regl,
      //   `https://picsum.photos/${Math.min(640, width)}/${Math.min(640, height)}/?random`, // `http://localhost/api/v1/acquisitions/4/Gd155/image?color=00ff5a&min=0&max=321`,
      //   true
      // );

      const img = new Image();
      img.crossOrigin = "";
      img.src = image;
      img.onload = () => {
        this.scatterplot = createScatterplot({
        regl,
        canvas,
        width,
        height,
        opacity: 1,
        pointSize: 2,
        pointSizeSelected: 1,
        pointOutlineWidth: 1,
        lassoMinDelay: 15,
        backgroundImage: regl.texture(img),
      });

      // const img = new Image();
      // img.crossOrigin = "";
      // img.src = "http://localhost/api/v1/acquisitions/4/Eu153/image?color=ff0000&min=0&max=159";
      // img.onload = () => {
      //   const regl = createRegl(canvas)
      //   this.scatterplot.set({ backgroundImage: regl.texture(img) })
      // };

      this.scatterplot.subscribe("pointover", this.pointoverHandler);
      this.scatterplot.subscribe("pointout", this.pointoutHandler);
      this.scatterplot.subscribe("select", this.selectHandler);
      this.scatterplot.subscribe("deselect", this.deselectHandler);

      window.addEventListener("resize", this.resizeHandler);

      // this.points = data.heatmap
      //   ? data.x.data.map((x, i) => {
      //       return [x, data.y.data[i], data.heatmap!.data[i], data.cell_ids[i]];
      //     })
      //   : data.x.data.map((x, i) => {
      //       return [x, data.y.data[i], 0, data.cell_ids[i]];
      //     });
      //
      // const categories = data.heatmap ? data.heatmap!.data : [0];
      // const min = Math.min(...categories);
      // const max = Math.max(...categories);
      // const colorScale = scaleSequential(interpolateCool).domain([min, max]);
      // const colors = [...Array(categories.length).keys()].map((item) => rgb(colorScale(item)).hex());
      // this.scatterplot.set({
      //   colorBy: "category",
      //   colors: colors,
      // });
      //
      // const p = this.points.map((item) => [item[0], item[1], item[2]]);
      // this.scatterplot.draw(p);

      const p = new Array(1000).fill(0).map(() => [
        -1 + Math.random() * 2, // x
        -1 + Math.random() * 2, // y
        "1_111", // category
        Math.random(), // value
      ]);
      this.scatterplot.draw(p);
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
        this.scatterplot.set({ width: rect.width, height: rect.height });
      }
    }
  }

  async mounted() {
    const canvas = this.$refs.canvas as any;

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
    this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);

    window.addEventListener("resize", this.resizeHandler);
  }

  beforeDestroy() {
    if (this.scatterplot) {
      this.scatterplot.destroy();
    }
  }
}
</script>
