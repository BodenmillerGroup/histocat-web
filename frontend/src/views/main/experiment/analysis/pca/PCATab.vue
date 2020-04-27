<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <v-chart :options="options" autoresize @brushselected="brushselected" />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>PCA Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
            <v-chip v-for="item in channels" :key="item" :value="item" small>
              {{ item }}
            </v-chip>
          </v-chip-group>
          <v-card-actions>
            <v-btn @click="selectAll" small :disabled="selectedChannels.length === channels.length">
              Select all
            </v-btn>
            <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0">
              Clear all
            </v-btn>
          </v-card-actions>
          <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
            <v-radio label="2D" value="2" />
            <v-radio label="3D" value="3" />
          </v-radio-group>
          <v-select
            :items="heatmaps"
            v-model="heatmap"
            label="Heatmap"
            hint="Heatmap marker"
            item-text="label"
            return-object
            persistent-hint
            clearable
            dense
            class="mt-5"
          />
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IPCAData } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import * as echarts from "echarts";
import "echarts-gl";
import "echarts/lib/chart/line";
import "echarts/lib/chart/scatter";
import "echarts/lib/component/toolbox";
import "echarts/lib/component/tooltip";
import "echarts/lib/component/visualMap";
import { uniq } from "rambda";
import { Component, Vue, Watch } from "vue-property-decorator";

const commonOptions: echarts.EChartOption = {
  title: {
    text: "Principal Component Analysis",
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
      brush: {
        title: ["Rect", "Polygon", "Clear"],
      },
    },
  },
  brush: {
    outOfBrush: {
      color: "#abc",
    },
    brushStyle: {
      borderWidth: 2,
      color: "rgba(0,0,0,0.2)",
      borderColor: "rgba(0,0,0,0.5)",
    },
    seriesIndex: [0, 1],
    throttleType: "debounce",
    throttleDelay: 1000,
    toolbox: ["rect", "polygon", "clear"],
  },
};

@Component
export default class PCATab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  options: echarts.EChartOption = {};

  selectedChannels: any[] = [];
  nComponents = "2";
  heatmap: { type: string; label: string } | null = null;

  points: any[] = [];

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    return this.showOptions ? 9 : 12;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get channels() {
    return this.datasetContext.getters.channels;
  }

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  selectAll() {
    this.selectedChannels = this.channels;
  }

  clearAll() {
    this.selectedChannels = [];
  }

  async submit() {
    let heatmap = "";
    if (this.heatmap) {
      heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
    }

    const acquisitionIds =
      this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

    await this.analysisContext.actions.getPCAData({
      dataset_id: this.activeDataset!.id,
      acquisition_ids: acquisitionIds,
      n_components: parseInt(this.nComponents, 10),
      heatmapType: this.heatmap ? this.heatmap.type : "",
      heatmap: heatmap,
      markers: this.selectedChannels,
    });
  }

  get pcaData() {
    return this.analysisContext.getters.pcaData;
  }

  @Watch("pcaData")
  pcaDataChanged(data: IPCAData) {
    if (data) {
      if (data.z) {
        this.plot3D(data);
      } else {
        this.plot2D(data);
      }
    }
  }

  private plot2D(data: IPCAData) {
    this.points = data.heatmap
      ? data.x.data.map((x, i) => {
          return [x, data.y.data[i], data.heatmap!.data[i]];
        })
      : data.x.data.map((x, i) => {
          return [x, data.y.data[i], data.cell_ids[i]];
        });

    const options: echarts.EChartOption = {
      ...commonOptions,
      xAxis: {
        type: "value",
        name: data.x.label,
        nameTextStyle: {
          fontWeight: "bold",
        },
      },
      yAxis: {
        type: "value",
        name: data.y.label,
        nameTextStyle: {
          fontWeight: "bold",
        },
      },
      dataset: {
        source: this.points,
        dimensions: [
          { name: data.x.label, type: "float" },
          { name: data.y.label, type: "float" },
        ],
      },
      series: [
        {
          type: "scatter",
          name: "Scatter2D",
          symbolSize: 2,
          large: false,
          encode: {
            x: 0,
            y: 1,
          },
          progressive: 0,
        },
      ],
    };

    if (data.heatmap) {
      (options.dataset as any).dimensions.push({ name: data.heatmap.label });
      options.visualMap = this.getVisualMap(data);
    }

    this.options = options;
  }

  private plot3D(data: IPCAData) {
    const options = {
      ...commonOptions,
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

  private getVisualMap(data: IPCAData): echarts.EChartOption.VisualMap[] {
    return data.heatmap!.label === "Acquisition"
      ? this.getCategoricalVisualMap(data)
      : this.getContinuousVisualMap(data);
  }

  private getCategoricalVisualMap(data: IPCAData): echarts.EChartOption.VisualMap[] {
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

  private getContinuousVisualMap(data: IPCAData): echarts.EChartOption.VisualMap[] {
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

  async brushselected(params) {
    const dataIndex = params.batch[0].selected[0].dataIndex;
    if (dataIndex.length > 0) {
      const cell_ids: number[] = [];
      for (const i of dataIndex) {
        const cell_id = Number(this.points[i][2].split("_")[1])
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
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
