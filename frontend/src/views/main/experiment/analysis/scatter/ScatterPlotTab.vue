<template>
  <v-layout row>
    <v-flex :class="mainClass">
      <v-chart
        :options="options"
        autoresize
        class="scatter-plot-view"
      />
    </v-flex>
    <v-flex
      v-if="showOptions"
      md3
    >
      <v-card tile>
        <v-card-title>Scatter Plot Settings</v-card-title>
        <v-card-text>
          <v-form
            v-model="valid"
            ref="form"
          >
            <v-select
              class="input-row"
              :items="items"
              item-text="label"
              v-model="markerX"
              return-object
              label="X"
              hint="X axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              class="input-row"
              :items="items"
              item-text="label"
              v-model="markerY"
              return-object
              label="Y"
              hint="Y axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              class="input-row"
              :items="items"
              item-text="label"
              v-model="markerZ"
              return-object
              label="Z"
              hint="Z axis marker"
              persistent-hint
            ></v-select>
            <v-switch
              v-model="showRegressionLine"
              label="Regression line"
              hide-details
            ></v-switch>
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-btn
            @click="submit"
            color="primary"
            block
            :disabled="!valid"
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
  import { IScatterPlotData } from '@/modules/analysis/models';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { IChannel } from '@/modules/experiment/models';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import { required } from '@/utils';
  import * as echarts from 'echarts';
  import ecStat from 'echarts-stat';
  import 'echarts/lib/chart/line';
  import 'echarts/lib/chart/scatter';
  import 'echarts/lib/component/toolbox';
  import 'echarts/lib/component/tooltip';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class ScatterPlotTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    readonly required = required;

    valid = false;

    options: echarts.EChartOption = {};

    showRegressionLine = false;
    markerX: IChannel | null = null;
    markerY: IChannel | null = null;
    markerZ?: IChannel | null = null;

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
      return this.activeAcquisition && this.activeAcquisition.channels ? this.activeAcquisition.channels : [];
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

        await this.analysisContext.actions.getScatterPlotData({
          datasetId: this.activeDataset.id,
          acquisitionId: this.activeAcquisition.id,
          markerX: this.markerX!.meta['OrderNumber'],
          markerY: this.markerY!.meta['OrderNumber'],
          markerZ: undefined,
        });
      }
    }

    get scatterPlotData() {
      return this.analysisContext.getters.scatterPlotData;
    }

    @Watch('scatterPlotData')
    scatterPlotDataChanged(data: IScatterPlotData) {
      if (!data) {
        return;
      }

      const points = data.x.map((x, i) => {
        return [x, data.y[i]];
      });

      const series: any = [
        {
          type: 'scatter',
          name: 'scatter',
          large: true,
          symbolSize: 4,
          encode: {
            // Map dimension "amount" to the X axis.
            x: 'X',
            // Map dimension "product" to the Y axis.
            y: 'Y',
            tooltip: ['X', 'Y'],
          },
        },
      ];

      if (this.showRegressionLine) {
        const myRegression = ecStat.regression('linear', points, 1);
        myRegression.points.sort((a, b) => {
          return a[0] - b[0];
        });
        series.push({
          name: 'line',
          type: 'line',
          large: true,
          showSymbol: false,
          data: myRegression.points,
        });
      }

      this.options = {
        animation: false,
        xAxis: {},
        yAxis: {},
        dataset: {
          source: points,
          dimensions: [
            { name: 'X', type: 'float' },
            { name: 'Y', type: 'float' },
          ],
        },
        series: series,
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
  .scatter-plot-view {
    height: calc(100vh - 212px);
  }

  .input-row {
    margin-bottom: 32px;
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
