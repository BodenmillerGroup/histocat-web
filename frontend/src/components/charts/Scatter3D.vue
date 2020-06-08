<template>
  <v-chart :options="options" autoresize class="chart" />
</template>

<script lang="ts">
import { IChart3DData } from "@/modules/analysis/models";
import * as echarts from "echarts";
import "echarts-gl";
import "echarts/lib/component/toolbox";
import "echarts/lib/component/tooltip";
import "echarts/lib/component/visualMap";
import { uniq } from "rambda";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";

@Component
export default class Scatter3D extends Vue {
  @Prop(Object) data;
  @Prop(String) title;

  options: echarts.EChartOption = {};

  @Watch("data")
  dataChanged(data: IChart3DData) {
    if (data) {
      const options = {
        title: {
          text: this.title,
          left: "center",
          top: 0,
        },
        animation: false,
        tooltip: {
          show: false,
        },
        toolbox: {
          show: true,
          right: "9%",
          feature: {
            restore: {
              show: true,
              title: "Reset",
            },
            saveAsImage: {
              show: true,
              title: "Export",
            },
            dataView: {
              show: true,
              title: "Data",
              readOnly: true,
              lang: ["Data View", "Hide", "Refresh"],
            },
          },
        },
        grid3D: {},
        xAxis3D: {
          type: "value",
          name: data.x.label,
          nameTextStyle: {
            fontWeight: "bold",
          },
        },
        yAxis3D: {
          type: "value",
          name: data.y.label,
          nameTextStyle: {
            fontWeight: "bold",
          },
        },
        zAxis3D: {
          type: "value",
          name: data.z!.label,
          nameTextStyle: {
            fontWeight: "bold",
          },
        },
        dataset: {
          source: [data.x.data, data.y.data, data.z!.data],
          dimensions: [
            { name: data.x.label, type: "float" },
            { name: data.y.label, type: "float" },
            { name: data.z!.label, type: "float" },
          ],
        },
        series: [
          {
            type: "scatter3D",
            name: "Scatter3D",
            seriesLayoutBy: "row",
            symbolSize: 2,
            encode: {
              x: data.x.label,
              y: data.y.label,
              z: data.z!.label,
              tooltip: [data.x.label, data.y.label, data.z!.label],
            },
          },
        ],
      } as echarts.EChartOption;

      if (data.heatmap) {
        (options.dataset as any).dimensions.push({ name: data.heatmap.label });
        (options.dataset as any).source.push(data.heatmap.data);
        options.visualMap = this.getVisualMap(data);
      }

      this.options = options;
    }
  }

  private getVisualMap(data: IChart3DData): echarts.EChartOption.VisualMap[] {
    return data.heatmap!.label === "Acquisition"
      ? this.getCategoricalVisualMap(data)
      : this.getContinuousVisualMap(data);
  }

  private getCategoricalVisualMap(data: IChart3DData): echarts.EChartOption.VisualMap[] {
    const categories = uniq(data.heatmap!.data);
    return [
      {
        type: "piecewise",
        orient: "vertical",
        top: "top",
        left: "right",
        categories: categories as any,
        padding: [
          60, // up
          20, // right
          5, // down
          5, // left
        ],
        inRange: {
          color: [
            "#e6194b",
            "#3cb44b",
            "#ffe119",
            "#4363d8",
            "#f58231",
            "#911eb4",
            "#46f0f0",
            "#f032e6",
            "#bcf60c",
            "#fabebe",
            "#008080",
            "#e6beff",
            "#9a6324",
            "#fffac8",
            "#800000",
            "#aaffc3",
            "#808000",
            "#ffd8b1",
            "#000075",
            "#808080",
            "#000000",
          ],
        },
      },
    ];
  }

  private getContinuousVisualMap(data: IChart3DData): echarts.EChartOption.VisualMap[] {
    const min = Math.min(...data.heatmap!.data);
    const max = Math.max(...data.heatmap!.data);
    return [
      {
        type: "continuous",
        orient: "horizontal",
        left: "center",
        text: ["Max", "Min"],
        calculable: true,
        realtime: false,
        min: min,
        max: max,
        inRange: {
          color: ["#4457cc", "#ff5200"],
        },
      },
    ];
  }
}
</script>

<style scoped>
.chart {
  height: calc(100vh - 88px);
  width: 100%;
}
</style>
