<template>
  <v-chart
    :options="options"
    autoresize
  />
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import * as echarts from 'echarts';
  import 'echarts/lib/chart/scatter';
  import 'echarts/lib/component/toolbox';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class ScatterPlotView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    options: echarts.EChartOption = {};

    get activeDataset() {
      return this.datasetContext.getters.activeDataset;
    }

    async mounted() {
      if (this.activeDataset) {
        const data = await this.analysisContext.actions.getScatterPlotData(this.activeDataset.id);
        this.options = {
          animation: false,
          xAxis: {},
          yAxis: {},
          dataset: {
            source: [
              data.x,
              data.y,
            ],
            dimensions: [
              { name: 'X', type: 'int' },
              { name: 'Y', type: 'float' },
            ],
          },
          series: [
            {
              type: 'scatter',
              name: 'scatter',
              seriesLayoutBy: 'row',
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
          ],
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
  }
</script>

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
