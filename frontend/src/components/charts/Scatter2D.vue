<template>
  <canvas ref="canvas" :width="width" :height="height" />
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import {
  interpolateCool,
  scaleSequential,
  scaleLinear,
  rgb,
  interpolateRainbow,
  interpolateSpectral,
  schemeCategory10, schemeTableau10
} from "d3";
import { experimentModule } from "@/modules/experiment";
import { selectionModule } from "@/modules/selection";
import { uniq } from "rambda";

@Component
export default class Scatter2D extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;
  @Prop(Number) width;
  @Prop(Number) height;

  points: any[] = [];
  scatterplot: any;
  selection: any[] = [];

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      this.scatterplot.deselect();
      // this.scatterplot.reset();

      this.points = data.heatmap
        ? data.x.data.map((x, i) => {
            return [x, data.y.data[i], data.heatmap!.data[i], data.acquisitionIds[i], data.cellIds[i]];
          })
        : data.x.data.map((x, i) => {
            return [x, data.y.data[i], data.acquisitionIds[i], data.acquisitionIds[i], data.cellIds[i]];
          });

      const categories = data.heatmap ? data.heatmap!.data : data.acquisitionIds;
      const min = Math.min(...categories);
      const max = Math.max(...categories);
      const colorScale = scaleSequential(interpolateRainbow).domain([min, max]);
      const colors: string[] = [];
      for (let i = 0; i < categories.length; i++) {
        colors.push(rgb(colorScale(i)).hex());
      }
      console.log(min, max);
      console.log(categories);
      console.log(uniq(colors));

      this.scatterplot.set({
        colorBy: "category",
        pointColor: schemeTableau10,
      });

      const p = this.points.map((item) => [item[0], item[1], item[2]]);
      this.scatterplot.draw(p);
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
      const cell_ids = new Map<number, number[]>();
      for (const i of this.selection) {
        const point = this.points[i];
        const acquisitionId = point[3];
        const cellId = point[4];
        if (!cell_ids.has(acquisitionId)) {
          cell_ids.set(acquisitionId, []);
        }
        const ids = cell_ids.get(acquisitionId);
        ids!.push(cellId);
      }
      this.selectionContext.mutations.setCellIds(cell_ids);
      if (this.applyMask) {
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
        this.scatterplot.set({ width: rect.width, height: rect.height });
      }
    }
  }

  async mounted() {
    const canvas = this.$refs.canvas as any;

    const { width, height } = canvas.getBoundingClientRect();

    // const regl = createRegl(canvas);
    // const backgroundImage = await createTextureFromUrl(
    //   regl,
    //   `https://picsum.photos/${Math.min(640, width)}/${Math.min(640, height)}/?random`, // `http://localhost/api/v1/acquisitions/4/Gd155/image?color=00ff5a&min=0&max=321`,
    //   true
    // );

    this.scatterplot = createScatterplot({
      // regl,
      canvas,
      width,
      height,
      opacity: 1,
      pointSize: 2,
      pointSizeSelected: 1,
      pointOutlineWidth: 1,
      lassoMinDelay: 15,
      // backgroundImage: backgroundImage,
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
  }

  beforeDestroy() {
    this.scatterplot.destroy();
  }
}
</script>
