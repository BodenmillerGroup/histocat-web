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
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class ScatterPlotView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    data = [];

    options: echarts.EChartOption = {
      xAxis: {},
      yAxis: {},
      series: [{
        symbolSize: 2,
        data: this.data,
        type: 'scatter',
      }],
    };

    get activeDataset() {
      return this.datasetContext.getters.activeDataset;
    }

    async mounted() {
      if (this.activeDataset) {
        const data = await this.analysisContext.actions.getScatterPlotData(this.activeDataset.id);
        const points = data.x.map((x, i) => {
          return [x, data.y[i]];
        });
        this.options = {
          xAxis: {},
          yAxis: {},
          series: [{
            symbolSize: 2,
            data: points,
            type: 'scatter',
          }],
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
