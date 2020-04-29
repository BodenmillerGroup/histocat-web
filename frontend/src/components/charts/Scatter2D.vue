<template>
  <canvas ref="canvas" :width="width" :height="height" />
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { interpolateCool, scaleSequential, rgb } from "d3";
import {experimentModule} from "@/modules/experiment";

@Component
export default class Scatter2D extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

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

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      this.scatterplot.deselect();
      this.scatterplot.reset();

      this.points = data.heatmap
        ? data.x.data.map((x, i) => {
            return [x, data.y.data[i], data.heatmap!.data[i], data.cell_ids[i]];
          })
        : data.x.data.map((x, i) => {
            return [x, data.y.data[i], 0, data.cell_ids[i]];
          });

      const categories = data.heatmap ? data.heatmap!.data : [0];
      const min = Math.min(...categories);
      const max = Math.max(...categories);
      const colorScale = scaleSequential(interpolateCool).domain([min, max]);
      const colors = [...Array(categories.length).keys()].map((item) => rgb(colorScale(item)).hex());
      this.scatterplot.set({
        colorBy: "category",
        colors: colors,
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
      const cell_ids: number[] = [];
      for (const i of this.selection) {
        const cell_id = Number(this.points[i][3].split("_")[1]);
        cell_ids.push(cell_id);
      }
      if (this.applyMask) {
        // this.settingsContext.mutations.setMaskSettings({
        //   ...this.settingsContext.getters.maskSettings,
        //   cell_ids: value,
        // });
        this.experimentContext.actions.getGatedMaskImage(cell_ids);
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

  mounted() {
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
    this.scatterplot.destroy();
  }
}
</script>
