<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <v-chart :options="options" autoresize />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>UMAP Settings</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
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
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  max="200"
                  step="1"
                  label="Neighbors"
                  v-model.number="nNeighbors"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="0.0"
                  max="0.99"
                  step="0.01"
                  label="Minimum distance"
                  v-model.number="minDist"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
            <v-select :items="metrics" v-model="metric" label="Metric" dense hide-details class="mt-5" />
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
        <v-card-text>
          <v-select
            :items="results"
            v-model="result"
            item-text="name"
            return-object
            label="Results"
            hint="UMAP processed data"
            persistent-hint
            clearable
            dense
            @change="resultChanged"
          />
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
          <v-btn @click="display" color="primary" block :disabled="!result">
            Display
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IUMAPData } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { required } from "@/utils/validators";
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
    text: "Uniform Manifold Approximation and Projection",
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
export default class UMAPTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;
  readonly metrics = [
    "euclidean",
    "manhattan",
    "chebyshev",
    "minkowski",
    "canberra",
    "braycurtis",
    "haversine",
    "mahalanobis",
    "wminkowski",
    "seuclidean",
    "cosine",
    "correlation"
  ];

  valid = false;

  options: echarts.EChartOption = {};

  selectedChannels: any[] = [];
  nComponents = "2";
  nNeighbors = 15;
  minDist = 0.1;
  metric = "euclidean";

  heatmap: { type: string; label: string } | null = null;

  result: any = null;

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

  get results() {
    return this.datasetContext.getters.umapResults;
  }

  selectAll() {
    this.selectedChannels = this.channels;
  }

  clearAll() {
    this.selectedChannels = [];
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const acquisitionIds =
        this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

      await this.analysisContext.actions.submitUMAP({
        dataset_id: this.activeDataset!.id,
        acquisition_ids: acquisitionIds,
        n_components: parseInt(this.nComponents, 10),
        markers: this.selectedChannels,
        n_neighbors: this.nNeighbors,
        min_dist: this.minDist,
        metric: this.metric
      });
    }
  }

  resultChanged(result) {
    if (result) {
      this.nComponents = result.params.n_components.toString();
      this.nNeighbors = result.params.n_neighbors;
      this.minDist = result.params.min_dist;
      this.metric = result.params.metric;
      this.selectedChannels = result.params.markers;

      this.experimentContext.mutations.setSelectedAcquisitionIds(result.params.acquisition_ids);
    }
  }

  async display() {
    if (!this.activeDataset) {
      self.alert("Please select a dataset");
      return;
    }

    if (!this.result) {
      self.alert("Please select result data");
      return;
    }

    let heatmap = "";
    if (this.heatmap) {
      heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
    }

    await this.analysisContext.actions.getUMAPResult({
      datasetId: this.activeDataset.id,
      name: this.result ? this.result.name : "",
      heatmapType: this.heatmap ? this.heatmap.type : "",
      heatmap: heatmap
    });
  }

  get umapData() {
    return this.analysisContext.getters.umapData;
  }

  @Watch("umapData")
  umapDataChanged(data: IUMAPData) {
    if (data) {
      if (data.z) {
        this.plot3D(data);
      } else {
        this.plot2D(data);
      }
    }
  }

  private plot2D(data: IUMAPData) {
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
        dimensions: [
          { name: data.x.label, type: "float" },
          { name: data.y.label, type: "float" }
        ]
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
      (options.dataset as any).dimensions.push({ name: data.heatmap.label });
      options.visualMap = this.getVisualMap(data);
    }

    this.options = options;
  }

  private plot3D(data: IUMAPData) {
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
      (options.dataset as any).dimensions.push({ name: data.heatmap.label });
      (options.dataset as any).source.push(data.heatmap.data);
      options.visualMap = this.getVisualMap(data);
    }

    this.options = options;
  }

  private getVisualMap(data: IUMAPData): echarts.EChartOption.VisualMap[] {
    return data.heatmap!.label === "Acquisition"
      ? this.getCategoricalVisualMap(data)
      : this.getContinuousVisualMap(data);
  }

  private getCategoricalVisualMap(data: IUMAPData): echarts.EChartOption.VisualMap[] {
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
          5 // left
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
            "#000000"
          ]
        }
      }
    ];
  }

  private getContinuousVisualMap(data: IUMAPData): echarts.EChartOption.VisualMap[] {
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
          color: ["#4457cc", "#ff5200"]
        }
      }
    ];
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
