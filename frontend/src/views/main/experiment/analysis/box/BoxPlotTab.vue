<template>
  <v-layout
    row
    class="chart-container"
  >
    <v-flex :class="mainClass">
      <v-chart
        :options="options"
        autoresize
      />
    </v-flex>
    <v-flex
      v-if="showOptions"
      md3
    >
      <v-card tile>
        <v-card-title>Box Plot Settings</v-card-title>
        <v-card-text>
          <v-chip-group
            v-model="selectedItems"
            multiple
            column
            active-class="primary--text"
          >
            <v-chip
              v-for="item in items"
              :key="item"
              :value="item"
              pill
            >
              {{ item }}
            </v-chip>
          </v-chip-group>
        </v-card-text>
        <v-card-actions>
          <v-btn
            @click="selectAll"
            small
            :disabled="selectedItems.length === items.length"
          >
            Select all
          </v-btn>
          <v-spacer></v-spacer>
          <v-btn
            @click="clearAll"
            small
            :disabled="selectedItems.length === 0"
          >
            Clear all
          </v-btn>
        </v-card-actions>
        <v-card-actions>
          <v-btn
            @click="submit"
            color="primary"
            block
            :disabled="selectedItems.length === 0"
          >
            Visualize
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-flex>
  </v-layout>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { IPlotSeries } from '@/modules/analysis/models';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import * as echarts from 'echarts';
  import 'echarts/extension/dataTool';
  import 'echarts/lib/chart/boxplot';
  import 'echarts/lib/chart/scatter';
  import 'echarts/lib/component/toolbox';
  import 'echarts/lib/component/tooltip';
  import { Component, Vue, Watch } from 'vue-property-decorator';

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

    get mainClass() {
      if (this.showOptions) {
        return 'md9';
      }
      return 'md12';
    }

    get activeAcquisition() {
      return this.experimentContext.getters.activeAcquisition;
    }

    get activeDataset() {
      return this.datasetContext.getters.activeDataset;
    }

    get items() {
      return this.activeDataset && this.activeDataset.artifacts['channel_map'] ? Object.keys(this.activeDataset.artifacts['channel_map']) : [];
    }

    selectAll() {
      this.selectedItems = this.items;
    }

    clearAll() {
      this.selectedItems = [];
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        if (!this.activeDataset) {
          self.alert('Please select a dataset');
          return;
        }

        if (!this.activeAcquisition) {
          self.alert('Please select an acquisition');
          return;
        }

        await this.analysisContext.actions.getBoxPlotData({
          datasetId: this.activeDataset.id,
          acquisitionId: this.activeAcquisition.id,
          markers: this.selectedItems,
        });
      }
    }

    get boxPlotData() {
      return this.analysisContext.getters.boxPlotData;
    }

    @Watch('boxPlotData')
    boxPlotDataChanged(plotData: IPlotSeries[]) {
      const data = (echarts as any).dataTool.prepareBoxplotData(plotData.map(item => item.data));
      data.axisData = plotData.map(item => item.marker);

      this.options = {
        animation: true,
        title: [
          {
            text: 'upper: Q3 + 1.5 * IQR \nlower: Q1 - 1.5 * IQR',
            borderColor: '#999',
            borderWidth: 1,
            textStyle: {
              fontSize: 14,
            },
            left: '12%',
            top: '10%',
          },
        ],
        xAxis: {
          type: 'category',
          data: data.axisData,
          splitLine: {
            show: false,
          },
          splitArea: {
            show: false,
          },
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        yAxis: {
          type: 'value',
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        series: [
          {
            type: 'boxplot',
            name: 'BoxPlot',
            data: data.boxData,
          },
          {
            name: 'outlier',
            type: 'scatter',
            large: true,
            data: data.outliers,
          },
        ],
        tooltip: {
          show: true,
        },
        toolbox: {
          show: true,
          right: '9%',
          feature: {
            restore: {
              show: true,
              title: 'Reset',
            },
            saveAsImage: {
              show: true,
              title: 'Export',
            },
          },
        },
      };
    }
  }
</script>

<style scoped>
  .chart-container {
    height: calc(100vh - 212px);
  }
</style>

<style>
  /**
   * The default size is 600px×400px, for responsive charts
   * you may need to set percentage values as follows (also
   * don't forget to provide a size for the container).
   */
  .echarts {
    width: 100%;
    height: 100%;
  }
</style>
