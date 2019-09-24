<template>
  <v-banner v-if="!activeDataset || !activeAcquisition" icon="mdi-alert-circle-outline">
    Please select acquisition and dataset
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <v-chart :options="options" autoresize />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>Box Plot Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedItems" multiple column active-class="primary--text">
            <v-chip v-for="item in items" :key="item" :value="item" small>
              {{ item }}
            </v-chip>
          </v-chip-group>
        </v-card-text>
        <v-card-actions>
          <v-btn @click="selectAll" small :disabled="selectedItems.length === items.length">
            Select all
          </v-btn>
          <v-btn @click="clearAll" small :disabled="selectedItems.length === 0">
            Clear all
          </v-btn>
        </v-card-actions>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedItems.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IPlotSeries } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import * as echarts from "echarts";
import "echarts/extension/dataTool";
import "echarts/lib/chart/boxplot";
import "echarts/lib/chart/scatter";
import "echarts/lib/component/toolbox";
import "echarts/lib/component/tooltip";
import { Component, Vue, Watch } from "vue-property-decorator";

@Component
export default class BoxPlotTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  options: echarts.EChartOption = {};

  selectedItems: any[] = [];

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

  selectAll() {
    this.selectedItems = this.items;
  }

  clearAll() {
    this.selectedItems = [];
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

    await this.analysisContext.actions.getBoxPlotData({
      datasetId: this.activeDataset.id,
      acquisitionId: this.activeAcquisition.id,
      markers: this.selectedItems
    });
  }

  get boxPlotData() {
    return this.analysisContext.getters.boxPlotData;
  }

  @Watch("boxPlotData")
  boxPlotDataChanged(data: IPlotSeries[]) {
    const boxplotData = (echarts as any).dataTool.prepareBoxplotData(data.map(item => item.data));
    boxplotData.axisData = data.map(item => item.label);

    this.options = {
      animation: true,
      title: [
        {
          text: "upper: Q3 + 1.5 * IQR \nlower: Q1 - 1.5 * IQR",
          borderColor: "#999",
          borderWidth: 1,
          textStyle: {
            fontSize: 14
          },
          left: "12%",
          top: "10%"
        }
      ],
      xAxis: {
        type: "category",
        data: boxplotData.axisData,
        splitLine: {
          show: false
        },
        splitArea: {
          show: false
        },
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      yAxis: {
        type: "value",
        nameTextStyle: {
          fontWeight: "bold"
        }
      },
      series: [
        {
          type: "boxplot",
          name: "BoxPlot",
          data: boxplotData.boxData
        },
        {
          name: "outlier",
          type: "scatter",
          large: true,
          data: boxplotData.outliers
        }
      ],
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
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
