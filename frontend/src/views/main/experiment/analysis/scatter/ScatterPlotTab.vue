<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition" icon="mdi-alert-circle-outline">
    Please select acquisition
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <v-chart :options="options" autoresize />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>Scatter Plot Settings</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
            <v-select
              :items="items"
              v-model="markerX"
              label="X"
              hint="X axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              :items="items"
              v-model="markerY"
              label="Y"
              hint="Y axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              :items="items"
              v-model="markerZ"
              label="Z"
              hint="Z axis marker"
              persistent-hint
              clearable
            ></v-select>
            <v-select
              :items="heatmaps"
              v-model="heatmap"
              label="Heatmap"
              hint="Heatmap marker"
              item-text="label"
              return-object
              persistent-hint
              clearable
              class="input-row"
            ></v-select>
            <v-switch v-if="!markerZ" v-model="showRegression" label="Show regression"></v-switch>
            <v-select
              v-if="!markerZ"
              class="input-row"
              :items="regressionTypes"
              v-model="regressionType"
              label="Regression type"
              hide-details
            ></v-select>
            <v-text-field
              v-if="!markerZ && regressionType === 'polynomial'"
              type="number"
              min="2"
              step="1"
              label="Polynomial order"
              v-model.number="polynomialOrder"
              :rules="[required]"
              hide-details
            ></v-text-field>
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="!valid">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IScatterPlotData } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { required } from "@/utils/validators";
import * as echarts from "echarts";
import "echarts-gl";
import ecStat from "echarts-stat";
import "echarts/lib/chart/line";
import "echarts/lib/chart/scatter";
import "echarts/lib/component/toolbox";
import "echarts/lib/component/tooltip";
import "echarts/lib/component/visualMap";
import { Component, Vue, Watch } from "vue-property-decorator";

type RegressionType = "linear" | "exponential" | "logarithmic" | "polynomial";

const commonOptions: echarts.EChartOption = {
  title: {
    text: "Mean intensity",
    left: "center",
    top: 0
  },
  animation: false,
  tooltip: {
    show: true
  },
  toolbox: {
    show: true,
    right: "9%",
    feature: {
      restore: {
        show: true,
        title: "Reset"
      },
      saveAsImage: {
        show: true,
        title: "Export"
      },
      dataView: {
        show: true,
        title: "Data",
        readOnly: true,
        lang: ["Data View", "Hide", "Refresh"]
      }
    }
  }
};

@Component
export default class ScatterPlotTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;
  readonly regressionTypes: RegressionType[] = ["linear", "polynomial"];

  valid = false;

  options: echarts.EChartOption = {};

  showRegression = false;
  regressionType: RegressionType = "linear";
  polynomialOrder = 2;

  markerX: string | null = null;
  markerY: string | null = null;
  markerZ: string | null = null;
  heatmap: { type: string; label: string } | null = null;

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    return this.showOptions ? 9 : 12;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get items() {
    return this.activeDataset && this.activeDataset.input["channel_map"]
      ? Object.keys(this.activeDataset.input["channel_map"])
      : [];
  }

  async submit() {
    if (!this.activeDataset) {
      self.alert("Please select a dataset");
      return;
    }

    if (!this.activeAcquisition) {
      self.alert("Please select an acquisition");
      return;
    }

    if ((this.$refs.form as any).validate()) {
      let heatmap = "";
      if (this.heatmap) {
        heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
      }

      await this.analysisContext.actions.getScatterPlotData({
        datasetId: this.activeDataset.id,
        acquisitionId: this.activeAcquisition.id,
        markerX: this.markerX!,
        markerY: this.markerY!,
        markerZ: this.markerZ ? this.markerZ : "",
        heatmapType: this.heatmap ? this.heatmap.type : "",
        heatmap: heatmap
      });
    }
  }

  get scatterPlotData() {
    return this.analysisContext.getters.scatterPlotData;
  }

  @Watch("scatterPlotData")
  scatterPlotDataChanged(data: IScatterPlotData) {
    if (data) {
      if (data.z) {
        this.plot3D(data);
      } else {
        this.plot2D(data);
      }
    }
  }

  private plot2D(data: IScatterPlotData) {
    const points = data.heatmap
      ? data.x.data.map((x, i) => {
          return [x, data.y.data[i], data.heatmap!.data[i]];
        })
      : data.x.data.map((x, i) => {
          return [x, data.y.data[i]];
        });

    const options: echarts.EChartOption = {
      ...commonOptions,
      xAxis: {
        type: "value",
        name: data.x.label,
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      yAxis: {
        type: "value",
        name: data.y.label,
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      dataset: {
        source: points,
        dimensions: [{ name: data.x.label, type: "float" }, { name: data.y.label, type: "float" }]
      },
      series: [
        {
          type: "scatter",
          name: "Scatter2D",
          symbolSize: 4,
          large: !data.heatmap,
          encode: {
            x: data.x.label,
            y: data.y.label,
            tooltip: [data.x.label, data.y.label]
          }
        }
      ]
    };

    if (data.heatmap) {
      (options.dataset as any).dimensions.push({ name: data.heatmap.label, type: "float" });

      options.visualMap = this.getVisualMap(data);
    }

    if (this.showRegression) {
      options.series!.push(this.getRegressionSeries(points));
    }

    this.options = options;
  }

  private plot3D(data: IScatterPlotData) {
    const options = {
      ...commonOptions,
      grid3D: {},
      xAxis3D: {
        type: "value",
        name: data.x.label,
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      yAxis3D: {
        type: "value",
        name: data.y.label,
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      zAxis3D: {
        type: "value",
        name: data.z!.label,
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      dataset: {
        source: [data.x.data, data.y.data, data.z!.data],
        dimensions: [
          { name: data.x.label, type: "float" },
          { name: data.y.label, type: "float" },
          { name: data.z!.label, type: "float" }
        ]
      },
      series: [
        {
          type: "scatter3D",
          name: "Scatter3D",
          seriesLayoutBy: "row",
          symbolSize: 4,
          encode: {
            x: data.x.label,
            y: data.y.label,
            z: data.z!.label,
            tooltip: [data.x.label, data.y.label, data.z!.label]
          }
        }
      ]
    } as echarts.EChartOption;

    if (data.heatmap) {
      (options.dataset as any).dimensions.push({ name: data.heatmap.label, type: "float" });

      (options.dataset as any).source.push(data.heatmap.data);

      options.visualMap = this.getVisualMap(data);
    }

    this.options = options;
  }

  private getVisualMap(data: IScatterPlotData): echarts.EChartOption.VisualMap[] {
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
          color: ["#4457cc", "#f45c00"]
        }
      }
    ];
  }

  private getRegressionSeries(points: ecStat.InputData) {
    const regression = ecStat.regression(this.regressionType, points, this.polynomialOrder);
    regression.points.sort((a, b) => {
      return a[0] - b[0];
    });
    return {
      name: "Regression",
      type: "line",
      lineStyle: {
        color: "#000000"
      },
      showSymbol: false,
      smooth: true,
      data: regression.points,
      markPoint: {
        itemStyle: {
          normal: {
            color: "transparent"
          }
        },
        label: {
          normal: {
            show: true,
            position: "left",
            formatter: regression.expression,
            textStyle: {
              color: "#000000",
              fontSize: 14
            }
          }
        },
        data: [
          {
            coord: regression.points[regression.points.length - 1]
          }
        ]
      }
    };
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}

.input-row {
  margin-bottom: 32px;
}
</style>
