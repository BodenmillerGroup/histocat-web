<template>
  <v-row no-gutters class="chart-container">
    <v-col
      :cols="columns"
    >
      <v-chart
        :options="options"
        autoresize
      />
    </v-col>
    <v-col
      v-if="showOptions"
      cols="3"
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
              v-model="markerX"
              label="X"
              hint="X axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              class="input-row"
              :items="items"
              v-model="markerY"
              label="Y"
              hint="Y axis marker"
              persistent-hint
              :rules="[required]"
            ></v-select>
            <v-select
              class="input-row"
              :items="items"
              v-model="markerZ"
              label="Z"
              hint="Z axis marker"
              persistent-hint
              clearable
            ></v-select>
            <v-select
              class="input-row"
              :items="heatmaps"
              v-model="heatmap"
              label="Heatmap"
              hint="Heatmap marker"
              item-text="label"
              item-value="value"
              persistent-hint
              clearable
            ></v-select>
            <v-switch
              v-if="!markerZ"
              v-model="showRegression"
              label="Show regression"
            ></v-switch>
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
    </v-col>
  </v-row>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { IScatterPlotData } from '@/modules/analysis/models';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import { required } from '@/utils';
  import * as echarts from 'echarts';
  import 'echarts-gl';
  import ecStat from 'echarts-stat';
  import 'echarts/lib/chart/line';
  import 'echarts/lib/chart/scatter';
  import 'echarts/lib/component/toolbox';
  import 'echarts/lib/component/tooltip';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  type RegressionType = 'linear' | 'exponential' | 'logarithmic' | 'polynomial';

  @Component
  export default class ScatterPlotTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    readonly required = required;
    readonly regressionTypes: RegressionType[] = ['linear', 'polynomial'];

    valid = false;

    options: echarts.EChartOption = {};

    showRegression = false;
    regressionType: RegressionType = 'linear';
    polynomialOrder = 2;

    markerX: string | null = null;
    markerY: string | null = null;
    markerZ: string | null = null;
    heatmap: string | null = null;

    get heatmaps() {
      return this.activeDataset && this.activeDataset.artifacts['neighbors_columns'] ? this.activeDataset.artifacts['neighbors_columns'].map(item => item.substring(10, item.length)) : [];
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
      return this.activeDataset && this.activeDataset.artifacts['channel_map'] ? Object.keys(this.activeDataset.artifacts['channel_map']) : [];
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
          markerX: this.markerX!,
          markerY: this.markerY!,
          markerZ: this.markerZ ? this.markerZ : '',
        });
      }
    }

    get scatterPlotData() {
      return this.analysisContext.getters.scatterPlotData;
    }

    @Watch('scatterPlotData')
    scatterPlotDataChanged(plotData: IScatterPlotData) {
      if (plotData) {
        if (plotData.z) {
          this.plot3D(plotData);
        } else {
          this.plot2D(plotData);
        }
      }
    }

    private plot2D(plotData: IScatterPlotData) {
      const points = plotData.x.data.map((x, i) => {
        return [x, plotData.y.data[i]];
      });

      const series: any = [
        {
          type: 'scatter',
          name: 'Scatter2D',
          large: true,
          symbolSize: 4,
          encode: {
            x: plotData.x.marker,
            y: plotData.y.marker,
            tooltip: [plotData.x.marker, plotData.y.marker],
          },
        },
      ];

      if (this.showRegression) {
        series.push(this.getRegressionSeries(points));
      }

      this.options = {
        animation: false,
        xAxis: {
          type: 'value',
          name: plotData.x.marker,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        yAxis: {
          type: 'value',
          name: plotData.y.marker,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        dataset: {
          source: points,
          dimensions: [
            { name: plotData.x.marker, type: 'float' },
            { name: plotData.y.marker, type: 'float' },
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

    private getRegressionSeries(points: ecStat.InputData) {
      const regression = ecStat.regression(this.regressionType, points, this.polynomialOrder);
      regression.points.sort((a, b) => {
        return a[0] - b[0];
      });
      return {
        name: 'Regression',
        type: 'line',
        showSymbol: false,
        smooth: true,
        data: regression.points,
        markPoint: {
          itemStyle: {
            normal: {
              color: 'transparent',
            },
          },
          label: {
            normal: {
              show: true,
              position: 'left',
              formatter: regression.expression,
              textStyle: {
                color: '#333',
                fontSize: 14,
              },
            },
          },
          data: [{
            coord: regression.points[regression.points.length - 1],
          }],
        },
      };
    }

    private plot3D(plotData: IScatterPlotData) {
      this.options = {
        animation: false,
        grid3D: {},
        xAxis3D: {
          type: 'value',
          name: plotData.x.marker,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        yAxis3D: {
          type: 'value',
          name: plotData.y.marker,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        zAxis3D: {
          type: 'value',
          name: plotData.z!.marker,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        dataset: {
          source: [
            plotData.x.data,
            plotData.y.data,
            plotData.z!.data,
          ],
          dimensions: [
            { name: plotData.x.marker, type: 'float' },
            { name: plotData.y.marker, type: 'float' },
            { name: plotData.z!.marker, type: 'float' },
          ],
        },
        series: [
          {
            type: 'scatter3D',
            name: 'Scatter3D',
            seriesLayoutBy: 'row',
            symbolSize: 4,
            encode: {
              x: plotData.x.marker,
              y: plotData.y.marker,
              z: plotData.z!.marker,
              tooltip: [plotData.x.marker, plotData.y.marker, plotData.z!.marker],
            },
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
      } as echarts.EChartOption;
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
