<template>
  <canvas ref="canvas" :width="canvasWidth" :height="canvasHeight" />
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
  schemeCategory10,
  schemeTableau10,
  interpolateReds,
} from "d3";
import { experimentModule } from "@/modules/experiment";
import { selectionModule } from "@/modules/selection";
import { forEach, uniq } from "rambda";
import { CellPoint } from "@/components/charts/CellPoint";

@Component
export default class Scatter2D extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;
  @Prop(Number) canvasWidth;
  @Prop(Number) canvasHeight;

  points: CellPoint[] = [];
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

      if (data.heatmap) {
        // Use heatmap value

        const values = data.heatmap!.data;
        const min = Math.min(...values);
        const max = Math.max(...values);

        const normalizedValues = values.map(v => v / max);

        this.points = data.x.data.map(
          (x, i) => new CellPoint(
            data.acquisitionIds[i],
            data.cellIds[i],
            x,
            data.y.data[i],
            normalizedValues[i])
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

  pointoverHandler(pointId: number) {
    const point = this.points[pointId];
    console.log(`X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nValue: ${point.value}`);
  }

  pointoutHandler(pointId: number) {
    const point = this.points[pointId];
    console.log(`X: ${point.x}\nY: ${point.y}\nAcquisitionId: ${point.acquisitionId}\nCellId: ${point.cellId}\nValue: ${point.value}`);
  }

  selectHandler({ points: selectedPoints }) {
    console.log("Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const cell_ids = new Map<number, number[]>();
      for (const i of this.selection) {
        const point = this.points[i];
        const acquisitionId = point.acquisitionId;
        if (!cell_ids.has(acquisitionId)) {
          cell_ids.set(acquisitionId, []);
        }
        const ids = cell_ids.get(acquisitionId);
        ids!.push(point.cellId);
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
    console.log(width, height)

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
    // this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);

    window.addEventListener("resize", this.resizeHandler);
  }

  beforeDestroy() {
    this.scatterplot.destroy();
  }
}
</script>
