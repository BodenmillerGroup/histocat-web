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
        <v-card-title>PhenoGraph Settings</v-card-title>
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
            <v-radio-group v-model="jaccard" mandatory hide-details label="Mode">
              <v-radio label="Jaccard metric" value="jaccard" />
              <v-radio label="Gaussian kernel" value="gaussian" />
            </v-radio-group>
            <v-row>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  max="200"
                  step="1"
                  label="Neighbors"
                  v-model.number="nearestNeighbors"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
              <v-col>
                <v-text-field
                  type="number"
                  min="2"
                  step="1"
                  label="Minimum cluster size"
                  v-model.number="minClusterSize"
                  :rules="[required]"
                  hide-details
                />
              </v-col>
            </v-row>
            <v-select :items="metrics" v-model="primaryMetric" label="Metric" hide-details dense class="mt-5" />
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
            hint="PhenoGraph processed data"
            persistent-hint
            clearable
            dense
            @change="resultChanged"
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
import { IPhenoGraphData } from "@/modules/analysis/models";
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
import { Component, Vue, Watch } from "vue-property-decorator";

@Component
export default class PhenoGraphTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;
  readonly metrics = ["euclidean", "manhattan", "correlation", "cosine"];

  valid = false;

  options: echarts.EChartOption = {};

  selectedChannels: any[] = [];
  nearestNeighbors = 30;
  jaccard = "jaccard";
  minClusterSize = 10;
  primaryMetric = "euclidean";

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

  get results() {
    return this.datasetContext.getters.phenographResults;
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

      await this.analysisContext.actions.submitPhenoGraph({
        dataset_id: this.activeDataset!.id,
        acquisition_ids: acquisitionIds,
        markers: this.selectedChannels,
        jaccard: this.jaccard === "jaccard",
        min_cluster_size: this.minClusterSize,
        nearest_neighbors: this.nearestNeighbors,
        primary_metric: this.primaryMetric
      });
    }
  }

  resultChanged(result) {
    if (result) {
      this.selectedChannels = result.params.markers;
      this.nearestNeighbors = result.params.nearest_neighbors;
      this.primaryMetric = result.params.primary_metric;
      this.minClusterSize = result.params.min_cluster_size;
      this.jaccard = result.params.jaccard ? "jaccard" : "gaussian";

      this.experimentContext.mutations.setSelectedAcquisitionIds(result.params.acquisition_ids);
    }
  }

  async display() {
    await this.analysisContext.actions.getPhenoGraphResult({
      datasetId: this.activeDataset!.id,
      name: this.result ? this.result.name : ""
    });
  }

  get phenographData() {
    return this.analysisContext.getters.phenographData;
  }

  @Watch("phenographData")
  phenographDataChanged(data: IPhenoGraphData) {
    if (data) {
      this.plot(data);
    }
  }

  private plot(data: any) {
    const communities = data.community;
    const markers = Object.keys(data).filter(item => item !== "community");
    const points: any[] = [];
    const mins: any[] = [];
    const maxs: any[] = [];
    for (let m = 0; m < markers.length; m++) {
      mins.push(Math.min(...data[markers[m]]));
      maxs.push(Math.max(...data[markers[m]]));
      for (let c = 0; c < communities.length; c++) {
        const v = data[markers[m]][c];
        points.push([c, m, v]);
      }
    }
    const options: echarts.EChartOption = {
      title: {
        text: "PhenoGraph",
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
      },
      xAxis: {
        type: "category",
        name: "Community",
        nameTextStyle: {
          fontWeight: "bold"
        },
        data: communities,
        splitArea: {
          show: true
        }
      },
      yAxis: {
        type: "category",
        name: "Marker",
        nameTextStyle: {
          fontWeight: "bold"
        },
        data: markers,
        splitArea: {
          show: true
        }
      },
      series: [
        {
          type: "heatmap",
          name: "PhenoGraph",
          data: points,
          itemStyle: {
            emphasis: {
              shadowBlur: 10,
              shadowColor: "rgba(0, 0, 0, 0.5)"
            }
          }
        }
      ],
      visualMap: [
        {
          min: Math.min(...mins),
          max: Math.max(...maxs),
          calculable: true,
          orient: "horizontal",
          left: "center"
        }
      ]
    };

    this.options = options;
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
